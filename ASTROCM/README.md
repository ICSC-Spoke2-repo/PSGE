### Shared pledge for the astro-particle community:
+ One objective of the PSGE project—besides the generalization of the [Computing Model (CM)](https://github.com/CYGNUS-RD/middleware/) through the CYGNO use case—is to standardize
access to pledged resources and optimize their utilization for small- and medium-scale experiments, thereby enabling efficient resource allocation 
and amortizing potential surges in computational demand.
+ To support this objective, an Identity and Access Management (IAM) service for authentication and authorization—currently limited to users with INFN platform credentials—has been
created and is managed by the CSN2 Computing Working Group ([GdL Calcolo](https://web.infn.it/csn2/index.php/it/struttura/gruppo-calcolo)): [iam.psge.cloud.infn.it](iam.psge.cloud.infn.it)
+ Through the [ui.psge.cloud.infn.it](ui.psge.cloud.infn.it) interface, notebook, VS Code, remote desktop, and remote SSH login environments are available, enabling users to develop, run, and manage workloads directly on PSGE resources.
+ The user interface includes a minimal software stack for analysis and access to CVMFS storage. Jobs can be submitted via the HTCondor CE endpoint tier2-ce2.lnf.infn.it; see the [PSGE Guide to Job Submission](./GuideJobSubmission.md)
 for step-by-step instructions.
