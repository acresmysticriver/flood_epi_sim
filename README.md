# Exploring the Conditional Quasi-Poisson Approach to Estimate Flooding Health Impacts with Simulated Data
This repository has been developed to demonstrate the process of creating a 'dummy' dataset with simulated event counts and exposure data in order to provide a system for testing modeling assumptions. The provided code is built to explore the impacts of flooding on health outcomes. Flooding impacts on health have been well established in the literature, but have often relied on single-event analysis. The growing availabilty of multimodal flood exposure data such as news catalogs, satellite products, stream gages, and meteorological observations provide the basis for more robust accounting of flood impacts across larger temporal and spatial extents. For a discussion of the available flood literature and an example approach to historical regional flood exposure assessment, see [Khemani et al.](https://iopscience.iop.org/article/10.1088/2752-5309/adedac/meta).

Code is separated into five sections that build iteratively toward a more complete pipeline for testing modeling assumptions. A summary of these sections is provided below:

## P1: Building a Dummy Event Count and Exposure Dataset
- Using the chicagoNMMAPS data structure embedded in the dlnm package, a daily event count time series is constructed with capacity for user-defined parameters of baseline count, variance, exposure effects, and time trend.
- Dummy flood events are inserted to the data structure randomly, with additional user-defined parameters of flood frequency and recurrence likelihood.
- The event and exposure data are then aggregated to a weekly event count and weekly marker of exposure.
