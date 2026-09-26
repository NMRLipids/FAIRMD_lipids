Simulation Curation
===================

Curation follows the workflow below. The first step is to determine whether the submission is
automated or manual: an automated pull request (PR) initially contains only ``info.yaml``, its name
starts with _"Upload Portal:"_ while a manual submission already contains the complete record,
including ``README.yaml`` metadata and computed properties (apl, eqtimes, OP, etc).

.. graphviz::

   digraph curation {
       graph [rankdir=TB, nodesep=0.35, ranksep=0.45, pad=0.1]
       node [shape=box, style="rounded", fontname="sans-serif", fontsize=16, margin="0.12,0.08"]
       edge [fontname="sans-serif", fontsize=13]

       A [label="Pull request"]
       B [label="Submission type?", shape=diamond]
       C [label="ValidateInfo check"]
       D [label="Validation successful?", shape=diamond]
       E [label="Fix metadata or reject PR"]
       F [label="Approve GitHub deployment"]
       G [label="Properties are computed"]
       H [label="Skip validation and deployment"]
       I [label="Check processing results"]
       J [label="Enrich metadata"]
       K [label="Problems?", shape=diamond]
       L [label="Fix and reprocess"]
       M [label="Close PR and request new upload"]
       N [label="Wait for checks and merge\ntemporary ID < 0"]
       O [label="UpdateID"]
       P [label="GlobalAnalysis"]
       Q [label="Occasional sanity check"]

       A -> B
       B -> C [label="Automated: info.yaml only"]
       C -> D
       D -> E [label="No"]
       E -> C
       D -> F [label="Yes"]
       F -> G
       B -> H [label="Manual: complete record"]
       G -> I
       H -> J
       I -> K
       K -> L [label="Fixable"]
       L -> I
       K -> M [label="New trajectory required"]
       K -> J [label="None"]
       J -> N
       N -> O -> P -> Q
   }

Automated submissions
---------------------

An automated submission starts when the PR contains only ``info.yaml``. The initial validation runs
with a 50 MB restriction, so it checks the metadata and processing setup without processing the full
trajectory. Typical problems include an incorrect mapping file, missing molecules in the composition, or
wrong residue names. The curator should fix these problems or ask the contributor to do so. If a
problem cannot be fixed in the repository, the PR should be rejected and the contributor should
correct the deposition before submitting again.

After validation, the curator manually approves the full processing GitHub deployment that will run
on the dedicated runner. Once approved, the PR indicates that it is being deployed. Processing
generates computed properties such as area per lipid, order parameters, form factor, and density
profiles, together with a nearly complete ``README.yaml``.

The curator then performs scientific sanity checks. In particular, the membrane should not have
collapsed, the equilibration time should be plausible, and the order parameters and other results
should not show obvious anomalies. Fixable problems in metadata should be corrected and checked
again. If a new trajectory upload is required, the current PR should be closed and the contributor
asked to create a new submission. Unfixable problems should be documented in a PR comment and the PR
rejected or closed; other contributors may be asked for assistance when appropriate.

The metadata bot will suggest corrections and additions, including Bioschema metadata. The curator
should accept these suggestions unless they contain obvious errors.

After all checks pass, merge the PR with a temporary ID (``< 0``). The post-merge ``UpdateID``
workflow assigns the next available unique ID, while ``GlobalAnalysis`` searches for
simulation--experiment matches, adds relationships, and calculates simulation quality where
applicable. These workflows do not require manual execution, but their results should be
sanity-checked occasionally.

Manual submissions
------------------

For a manual submission, processing and post-processing have already been performed by the
contributor, and the PR contains raw metadata together with computed properties. The curator skips
the ``info.yaml`` validation and deployment stages, then follows the remaining workflow: check the
record, enrich its metadata, wait for all checks, merge with a temporary ID, and let ``UpdateID`` and
``GlobalAnalysis`` complete the post-merge processing.

Troubleshooting
---------------

1. ``ValidateInfo`` passed but the actual processing failed. What to do?

    The curator should check the processing logs and the PR comments for any errors. If the problem is
    fixable, the curator can correct the metadata or processing setup and re-run the processing (by closing
    and re-opening the PR). If the problem is on the project's side, the curator should **open an issue** in
    the corresponding repository and **convert the PR to a draft** until the issue is resolved. The curator should also **comment on the PR** to inform the contributor about the problem and the next steps.

2. Issue is fixed, how to retrigger the process?

    If there are changes in workflows, the branch may need to be rebased. Please, synchronize ``UserData`` main
    branch to the ``BilayerData`` main and rebase the PR branch on top of it. Then, close and re-open the PR to retrigger the processing.

3. User made a new upload, can I update an old PR?

    If user made new upload in Zenodo, the curator can fix zenodo's DOI in the ``info.yaml`` and
    re-run the processing. If the user made a new upload using UploadPortal, it's a duplicate - one
    of them should be closed.
