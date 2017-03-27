#  The Change Log
This is a bit overdue, but better better Nate than lever, right?
----
## 0.3.0
### Features
- added filters to halt subassembly of contigs with poor coverage
- added CLI option for min_cov_depth, where pseudocontigs failing to reach a minimum of 5x coverage (by default) will e rejected from future iterations

### Fixes
- improved printing diagnostics

## 0.2.0
### Features
- changes to default parameters
### Fixes
- improved pipeline congruency
- fixed bug causing skipped iterations

## 0.1.0
This is the first logged change.  We changed the commandline arguments so that it now supports a single-end library as input, not just paired end data.
### Features
- supports single-end library
### Fixes
- fixed library parsing realting to setting `libtype`
