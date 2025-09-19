# Exploring the Conditional Quasi-Poisson Approach to Estimate Flooding Health Impacts with Simulated Data
This repository has been developed to demonstrate the process of creating a 'dummy' dataset with simulated event counts and exposure data in order to provide a system for testing modeling assumptions. The provided code is built to explore the impacts of flooding on health outcomes. Flooding impacts on health have been well established in the literature, but have often relied on single-event analysis. The growing availabilty of multimodal flood exposure data such as news catalogs, satellite products, stream gages, and meteorological observations provide the basis for more robust accounting of flood impacts across larger temporal and spatial extents. For a discussion of the available flood literature and an example approach to historical regional flood exposure assessment, see [Khemani et al.](https://iopscience.iop.org/article/10.1088/2752-5309/adedac/meta).

Code is separated into five sections that build iteratively toward a more complete pipeline for testing modeling assumptions. A summary of these sections is provided below:

## P1: Building a Dummy Event Count and Exposure Dataset
- Using the chicagoNMMAPS data structure embedded in the dlnm package, a daily event count time series is constructed with capacity for user-defined parameters of baseline count, variance, exposure effects, and time trend.
- Dummy flood events are inserted to the data structure randomly, with additional user-defined parameters of flood frequency and recurrence likelihood.
- The event and exposure data are then aggregated to a weekly event count and weekly marker of exposure.

## P2: Expanding Dummy Data and Inserting Flood Effects
- The data development process is transitioned to a function to allow for more seamless iteration through parameters and developing different inputs across dummy counties.
- For each county, user provided RR estimated for flooding are implemented to the event count time series to approximate flood impacts.
- The resulting event time series with exposure impacts is visualized across the 2 fake counties.

## P3: Adding Lag Effects and Selecting Controls for Conditional Quasi-Poisson Approach
- Script 1 is further modified to integrate a lag effect for flooding, where all four weeks following a flood also see an increase in event counts.
- Script 2 is added to select control weeks for comparison with the flood weeks defined in script 1.
- Following the approach of [Aggarwal et al.](https://arxiv.org/abs/2309.13142), control weeks are selected for the same week-of-year in the years preceding and following a flood event.
- The control selection process is iterative so as to only select control periods of 5 weeks (flood plus lag period) where no flooding was present.
