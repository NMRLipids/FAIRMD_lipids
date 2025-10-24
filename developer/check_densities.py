import json

import yaml

from fairmd.lipids.core import initialize_databank

if __name__ == "__main__":
    systems = initialize_databank()

    for system in systems:
        WaterDensity_name = system["path"] + "WaterDensity.json"
        try:
            f = open(WaterDensity_name)
            # print('Density file not found')
            WaterDensity = json.load(f)
            center = round(len(WaterDensity) / 2)
            # print(WaterDensity[center][1], WaterDensity[0][1])
            if WaterDensity[center][1] > WaterDensity[0][1]:
                print(system["path"])
                system["WARNINGS"] = {}
                system["WARNINGS"]["PBC"] = "z-jumps"
                # print(system)

                outfileDICT = system["path"] + "/README.yaml"
                with open(outfileDICT, "w") as f:
                    yaml.dump(system, f, sort_keys=False)

                # os.system('rm ' + system['path'] + '/*Density*')
                # os.system('rm ' + system['path'] + '/FormFactor*')
        except Exception:
            pass
            # print('Density file not found from ' +  system['path'])
