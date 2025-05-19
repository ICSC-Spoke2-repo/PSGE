# Guide on how to submit a Job to the PSGE Queue

---

**List of Available CEs at FRASCATI-IT (LNF):**

| Hostname                                    | vCPU | RAM (GiB) | OS          | Authorized VOs                                         |
| ------------------------------------------- | ---- | --------- | ----------- | ------------------------------------------------------ |
| **tier2-ce2.lnf.infn.it**                   |  400 | 1200      | Almalinux 9 | dteam, belle, ops, atlas, padme, lhcb, cta, csn2(psge) |

---

**Submission via Scitoken:**

Access the CSN2/PSGE jupyter notebook: [https://ui.psge.cloud.infn.it/](https://ui.psge.cloud.infn.it/), and configure the token:

```bash
[abrittac@ui ~]$ eval `oidc-agent`
[abrittac@ui ~]$ oidc-gen -w device
```

---

**IAM Provider Configuration:**

Select the issuer, e.g.:

* `https://iam.psge.cloud.infn.it/`

**Select the required scopes:**

```
openid compute.create offline_access profile compute.read compute.cancel compute.modify wlcg wlcg.groups
```

---

**Client Registration:**

Open the indicated link and enter the displayed code:

* [https://iam.psge.cloud.infn.it/device](https://iam.psge.cloud.infn.it/device)

---

**Completing the Configuration:**

```bash
[abrittac@ui ~]$ export BEARER_TOKEN=$(oidc-token test_lnf)
[abrittac@ui ~]$ mask=$(umask); umask 0077 ; oidc-token test_lnf > ${HOME}/token ; umask $mask
[abrittac@ui ~]$ export _condor_SEC_CLIENT_AUTHENTICATION_METHODS=SCITOKENS
```

---

**Check Write Permissions:**

```bash
condor_ping -verbose -pool tier2-ce2.lnf.infn.it:9619 -name tier2-ce2.lnf.infn.it -type schedd write
```

**Expected Output:** `Authorized: TRUE`

---

**Submitting a Test Job**

Create `test.sub` and `test.sh` files:

**test.sub**

```bash
universe = vanilla
scitokens_file = $ENV(HOME)/token
+owner = undefined
executable = test.sh
output = test.out
error = test.err
log = test.log
should_transfer_files = yes
when_to_transfer_output = on_exit
queue
```

**test.sh**

```bash
#!/bin/bash
hostname
whoami
pwd
cat /etc/redhat-release
```

---

**Submit to Tier2 LNF:**

```bash
[abrittac@ui ~]$ condor_submit --spool --name tier2-ce2.lnf.infn.it --pool tier2-ce2.lnf.infn.it:9619 test.sub
```

---

**Check Job Status**

```bash
[abrittac@ui ~]$ condor_q -pool tier2-ce2.lnf.infn.it:9619 -name tier2-ce2.lnf.infn.it <JOB ID>
```

---

**Retrieve Job Output**

```bash
[abrittac@ui ~]$ condor_transfer_data -pool tier2-ce2.lnf.infn.it:9619 -name tier2-ce2.lnf.infn.it <JOB ID>
```

---
**Enviroment variable that can make life easier (optional):**
```
CE=tier2-ce2.lnf.infn.it
export _condor_SEC_CLIENT_AUTHENTICATION_METHODS=SCITOKENS

export _condor_CONDOR_HOST=$CE:9619
export _condor_SCHEDD_HOST=$CE

alias condor_submit='condor_submit -spool'
```
