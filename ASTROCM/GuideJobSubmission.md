# Guide on how to submit a Job to the PSGE Queue

In PSGE seemsless job submission is enable by means of a centralized
HTCondor cluster. This cluster is currently powered by workers made
available by the Frascati (INFN-LNF) and Catania (INFN-CT) datacentres.

A set of INDIGO-IAM instances is integrated with the PSGE cluster,
allowing their users to exploit the native HTCondor token renewal via
credmon. These tokens, when requested and authorized by the user, are
made available and renewed through the job lifetime automatically and
can be used to access remote storage elements. See below for further
information.

## User access to PSGE

The job submission node is the PSGE user interface, that is accessible
via web through the jupyter hub interface,
[https://ui.psge.cloud.infn.it/](https://ui.psge.cloud.infn.it/), and
via SSH at the same address:

```bash
ssh <psge username>@ui.psge.cloud.infn.it
```

User access is always mediated via the PSGE identity and access
management system, i.e. [IAM-PSGE](https://iam.psge.cloud.infn.it/).

To enable ssh access, you need to upload your SSH public key in your
IAM-PSGE profile page.

## Job submission

First of all make sure to be able to open a terminal on the PSGE user
interface, either via SSH or via the terminal embedded in a jupyter
notebook.

Then, prepare an HTCondor submit file. Example can be found in the
[Official HTcondor user
guide](https://htcondor.readthedocs.io/en/latest/users-manual/submitting-a-job.html)
and in the [INFN-Tier 1 user
guide)[https://confluence.infn.it/spaces/TD/pages/75435014/Examples].

Here below an example is also reported, submitting **three times** (see
the last line) a simple shell script file, `script.sh`, that requires
**two cores** and enabling the automatic renewal of iam-cygno access
tokens (see second- and third-last lines) that can be used to access a
storage element from a running job.

```
executable              = script.sh
arguments               = 10m arg2 arg3 arg4
output                  = stdout-$(ClusterId).$(ProcID).txt
error                   = stderr-$(ClusterId).$(ProcID).txt
log                     = output-$(ClusterId).$(ProcID).log
request_cpus            = 2
transfer_executable     = Yes
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
want_io_proxy           = true
use_oauth_services      = cygno
cygno_oauth_permissions = profile,email,openid,offline_access,wlcg.groups
queue 3
```

The `script.sh` file being, for example, a simple hello world:

```bash
#!/bin/bash

echo 'hello world!'
echo "Arguments: $@"

sleep $1
```

Then, submit the jobs using the `condor_submit` command:

```bash
condor_submit -spool submit.sub
```

When you request the automatic renewal of access tokens, you may need to
execute the submit command twice. The first time may be to delegate
HTCondor the access token renewal. The delegation procedure is performed
by clicking on the link shown in the `condor_submit` output, for example:

```
[cpellegr@ui credmon]$ condor_submit -spool submit.sub
Submitting job(s)
Hellow, cpellegr.
Please visit: https://condor.psge.cloud.infn.it/key/XXXXXXXXXXXXXXREDACTEDXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```

This link is univocally linked to your identity. By clicking on it, a
login page is presented. Log into the chosen IAM instance, in the
example iam-cygno, and authorize the use of the listed scopes.

By successfully performing the login and scope authorization procedure,
you are delegating HTCondor to continuously provide your running jobs
with valid access tokens. This delegation lasts up to some days after
the last job ends and cannot excede 30 days.

By repeating the same job submission command, after having successfully
delegated HTCondor, you actually submit the jobs, for example:

```
[cpellegr@ui credmon]$ condor_submit -spool submit.sub
Submitting job(s)
Hellow, cpellegr.
Please visit: https://condor.psge.cloud.infn.it/key/XXXXXXXXXXXXXXREDACTEDXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

[cpellegr@ui credmon]$ condor_submit -spool submit.sub
Submitting job(s)...
3 job(s) submitted to cluster 43.
```

If you don't need to delegate HTCondor to renewing access tokens for
you, remove the

```
use_oauth_services      = cygno
cygno_oauth_permissions = profile,email,openid,offline_access,wlcg.groups
```

lines from the submit file.

These examples are for the iam-cygno instance, but the PSGE HTCondor can be and is enabled with different token issuers.

Currently the following IAM instances are enabled:

| IAM-Instance                        | name of the oauth service in PSGE | Description |
| ----------------------------------- | --------------------------------- | ----------- |
| iam-t1-computing.cloud.cnaf.infn.it | t1                                | multi-VO    |
| iam-cygno.cloud.cnaf.infn.it        | cygno                             | CYGNO-VO    |
| iam.cloud.infn.it                   | infncloud                         | multi-VO    |

To enable the `t1` one, for example, you need to modify both the
`use_oauth_services` line and the next one:

```
use_oauth_services   = t1
t1_oauth_permissions = profile,email,openid,offline_access
```
By replacing `cygno` with `t1` in both.

Enabling more than one delegation is possible, by adding more keys to
the `use_oauth_services` submit file line and one `*_oauth_permissions`
line per key, for example:

```
use_oauth_services      = t1,cygno
t1_oauth_permissions    = profile,email,openid,offline_access
cygno_oauth_permissions = profile,email,openid,offline_access,wlcg.groups
```

### Using the access token

To use the access token automatically renewed by the PSGE HTCondor you
have delegated from a job, you need to read the corresponding file,
always located in the `${_CONDOR_CREDS}` folder present in the job
environment. The files are in `JSON` format and the name is in the form
`<oauth_service_name>.use`. For example, if you used the `cygno` IAM
instance, the file will be located in `${_CONDOR_CREDS}/cygno.use` on
the worker node. Please note that this is an absolute file path, so it
remains valid despite any `cd` or `pushd` commands in the job script.

To load the access token as the classical `BEARER_TOKEN` environment
variable to be used, for example, with the gfal toolkit, you can use the
following command line:

```bash
export BEARER_TOKEN="$(jq -r .access_token "${_CONDOR_CREDS}/cygno.use")"
```

Make sure to always load a fresh token before attempting to access a
grid storage element:

```bash
export BEARER_TOKEN="$(jq -r .access_token "${_CONDOR_CREDS}/t1.use")"
gfal-copy my-huge-output.dat https://xfer-archive.cr.cnaf.infn.it:8443/swgo/...
```

