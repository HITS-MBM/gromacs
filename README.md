# GROMACS developments by Molecular Biomechanics (MBM) group at HITS

## Warning: development code

The GROMACS features available in a number of branches in this repository are a work in progress code being developed at HITS. The features are considered generally usable, being used in our group (and some of the collaborating groups) on a more or less regular basis.

However, the features as provided here should not be considered stable. They are prone to changes as needed i.e. we don't provide any compatibility guarantees for the future: the options can be added, removed, or extended as needed and the option names can be changed. Any of those changes can occur without the documentation describing the migration process from the previous version.

Furthermore, the feature branches can *and likely will* be rebased on top of latest upstream release/master versions as needed. To be safe, please use

```
git pull --rebase
```

when updating. The changes will be [submitted to GROMACS](https://gerrit.gromacs.org/) for the review and a possible inclusion when considered ready.

tl;dr **Use the code at your own risk.**

## Feature suggestions and bug reports

We welcome usage reports, suggestions for changes, and especially bug reports. Please provide the input files and the sequence of commands required to reproduce the issue.

## Features available

Features are provided in separate branches named according to the GROMACS version upon which they are based. Each of the features is briefly described in the `README.md` file of its respective branch.

### Conditional Stop

See branch `conditional-stop-release-2016`.

### Sliced Pull

See branches `sliced-pull-release-2016` and `sliced-pull-release-2018`.
