import os
import sys
import requests
import yaml

def check_api(url):
    try:
        resp = requests.get(url, timeout=5)
        return resp.status_code == 200
    except requests.RequestException:
        return False

def get_chembl(inchikey):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?standard_inchi_key={inchikey}&format=json"
    if check_api(url):
        resp = requests.get(url)
        return resp.json() if resp.status_code == 200 else {}
    return {}

def get_pubchem(inchikey):
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        f"{inchikey}/property/IUPACName,SMILES,InChI,InChIKey,MolecularFormula,MolecularWeight/JSON"
    )
    if check_api(url):
        resp = requests.get(url)
        if resp.status_code == 200:
            return resp.json()["PropertyTable"]["Properties"][0]
    return {}

def get_pubchem_synonyms(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    if check_api(url):
        resp = requests.get(url)
        if resp.status_code == 200:
            return resp.json().get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
    return []

def get_chebi(chebi_id):
    if not chebi_id:
        return {}
    
    url = f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{chebi_id}/?only_ontology_parents=false&only_ontology_children=false"
    
    try:
        if check_api(url):
            resp = requests.get(url)
            if resp.status_code == 200:
                return resp.json()
    except requests.RequestException:
        pass
    
    return {}

def get_unichem(inchikey):
    url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
    if check_api("https://www.ebi.ac.uk/unichem/api/v1/sources"):
        headers = {'Content-Type': 'application/json'}
        data = {"type": "inchikey", "compound": inchikey}
        resp = requests.post(url, headers=headers, json=data)
        if resp.status_code == 200:
            compounds = resp.json().get("compounds", [])
            if compounds and "sources" in compounds[0]:
                return compounds[0]["sources"]
    return []

def extract_sameas(sources):
    mapping = {
        "pubchem": "pubchem.compound",
        "chebi": "ChEBI",
        "chembl": "ChEMBL",
        "lipidmaps": "lipidmaps",
        "metabolights": "metabolights",
        "swisslipids": "slm",
        "pdb": "pdb.ligand",
        "unii": "unii",
        "cas": "cas"
    }
    result = {}
    for src in sources:
        prefix = mapping.get(src["shortName"])
        if prefix:
            value = src["compoundId"]
            if prefix == "ChEBI":
                value = f"CHEBI:{value}" if value else ""
            elif prefix == "pubchem.compound":
                try:
                    value = int(value)
                except ValueError:
                    pass
            result[prefix] = value
    return result

def get_chembl_id_from_unichem(sources):
    for src in sources:
        if src["shortName"] == "chembl":
            return src["compoundId"]
    return None

def load_existing_metadata(path):
    if os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    return {}

def update_metadata(existing, new_data):
    for key, value in new_data.items():
        if isinstance(value, dict):
            updated = update_metadata(existing.get(key, {}), value)
            if updated:
                existing[key] = updated
        elif isinstance(value, list):
            if value:
                existing[key] = existing.get(key, []) or value
        else:
            if value not in [None, "", {}]:
                existing[key] = existing.get(key) or value
    return existing

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <metadata.yaml path>")
        sys.exit(1)

    metadata_path = sys.argv[1]

    # Extract NMRlipidsID from path (assumes structure: Molecules/membrane/<NMRlipidsID>/metadata.yaml)
    try:
        nmr_id = os.path.basename(os.path.dirname(metadata_path))
    except Exception:
        print("Error: Could not extract NMRlipidsID from path.")
        sys.exit(1)

    existing = load_existing_metadata(metadata_path)
    try:
        inchikey = existing["bioschema_properties"]["inChIKey"]
    except Exception:
        print("Error: Could not find bioschema_properties -> inChIKey in YAML file.")
        sys.exit(1)

    chembl = get_chembl(inchikey)
    pubchem = get_pubchem(inchikey)
    sources = get_unichem(inchikey)
    sameas = extract_sameas(sources)

    cid = pubchem.get("CID", sameas.get("pubchem.compound"))
    synonyms = get_pubchem_synonyms(cid) if cid else []

    # First, check if there's a ChEBI ID from unichem
    chebi_id = sameas.get("ChEBI", "").replace("CHEBI:", "")
    chebi_data = get_chebi(chebi_id) if chebi_id else {}
    
    # Collect alternate names with priority
    alternate_names = []
    
    # 1. Try ChEBI synonyms first
    if chebi_data and 'names' in chebi_data:
        # Extract only the 'name' from SYNONYM type
        alternate_names = [
            syn['name'] 
            for syn in chebi_data.get('names', {}).get('SYNONYM', []) 
            if syn.get('type') == 'SYNONYM' and syn.get('name')
        ]
    
    # 2. If no ChEBI synonyms, try ChEMBL synonyms
    if not alternate_names and chembl.get("molecule_synonyms"):
        alternate_names = [
            syn.get('molecule_synonym', '') 
            for syn in chembl.get("molecule_synonyms", [])
        ]
    
    # 3. If still no synonyms, try PubChem synonyms
    if not alternate_names and synonyms:
        alternate_names = synonyms

    molecule_props = chembl.get("molecule_properties", {})
    molecule_structures = chembl.get("molecule_structures", {})
    chembl_id = get_chembl_id_from_unichem(sources)

    # Image selection logic
    if chembl_id:
        image_url = f"https://www.ebi.ac.uk/chembl/api/data/image/{chembl_id}?dimensions=200"
    elif cid:
        image_url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l"
    else:
        image_url = ""

    bioschema = {
        "name": molecule_props.get("iupac_name") or pubchem.get("IUPACName", ""),
        "iupacName": molecule_props.get("iupac_name") or pubchem.get("IUPACName", ""),
        "molecularFormula": molecule_props.get("full_molformula") or pubchem.get("MolecularFormula", ""),
        "molecularWeight": float(molecule_props.get("full_mwt") or pubchem.get("MolecularWeight", 0)),
        "inChI": molecule_structures.get("standard_inchi") or pubchem.get("InChI", ""),
        "inChIKey": molecule_structures.get("standard_inchi_key") or pubchem.get("InChIKey", ""),
        "smiles": molecule_structures.get("canonical_smiles") or pubchem.get("SMILES", ""),
        "image": image_url,
        "description": ""
    }

    if alternate_names:
        bioschema["alternateName"] = alternate_names

    new_data = {
        "NMRlipids": {
            "id": nmr_id,
            "name": "",
            "charge": ""
        },
        "sameAs": sameas,
        "bioschema_properties": bioschema
    }

    updated = update_metadata(existing, new_data)

    with open(metadata_path, 'w', encoding='utf-8') as f:
        yaml.dump(updated, f, sort_keys=False, allow_unicode=True, default_flow_style=False)

    print(f"Updated metadata written to {metadata_path}")

if __name__ == "__main__":
    main()
