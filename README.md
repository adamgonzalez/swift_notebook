# swift_notebook
A collection of `R` scripts to perform multi-wavelength light curve analysis from the *Swift* X-ray observatory. This is still a work in progress, as there are a lot of methods to incorporate!

The methods are all contained within `functions.R`, while the workbook walkthrough is contained in `notebook.md`.

To do:
- Incorporate the Savitzky-Golay de-trending
	- Re-run the SF and ICCF analyses with de-trended light curves
	- Explore the differences with previous results
- Perform the flux-flux analysis
	- Extract the variable and constant components
