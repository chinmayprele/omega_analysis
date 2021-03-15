## Gene Information

## Primary Authors

- J. Yordy


### Idiosyncrasies of Run

- Each windowed run should take about 5-8 seconds each with a few `time.sleep()` statements.
- However, a few runs ran for 130s+, which was unexpected.
	- Windows 6, 7, and 8 had this issue, but due to lack of constant monitoring, addition of other windows to this list is currently impossible.
	- There are certain runs that have very high omegas (and some with `NaN` omegas).
	- This might be the cause of the issue.
	- Furthermore, it could also be because "omegas (dN/dS)" could not be `grep`ped out of the file. This seems unlikely because near the `NaN`, the omegas are rather high, and the omegas could be artificially high. More analysis at this region needs to be incorporating more species.
