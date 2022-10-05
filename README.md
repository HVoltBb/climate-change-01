![cover image](https://hvoltbb.github.io/pics/blueshark.jpg)

# Tracking the effects of climate change on marine life - case study 1

![version: v0.5](https://img.shields.io/badge/version-v0.5-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Platform: win | linux | ios](https://img.shields.io/badge/platform-win%20%7C%20linux%20%7C%20ios-lightgrey)

## Quick View

**Background**: Substantial progress has been made in identifying large-scale climate effect on somatic growth through the use of ageing-based methods in aquatic environments, yet their annual/seasonal temporal resolution seems too coarse for such a fast process.

**Problem statement**: Temporal resolution is a missing dimension in our understanding of climate effects on growth.

**Contributions**: An alternative source of high temporal resolution growth increments embedded within a multi-decadal traditional tag-recapture database was analysed to identify climate signals in the somatic growth of blue sharks Prionace glauca in the North Atlantic. Results indicate the growth response of P. glauca to the NAO occurred at a daily scale with a time-lag. Non-parametric modelling reveals an optimal response curve around the historical average of the NAO, and a significant negative response for large positive NAO anomalies. The temporal resolution of this study is unprecedented among current ageing-based studies with a comparable temporal coverage. Integrating high temporal resolution into long-term climate effect studies can open new avenues for research on identifying climate effect on growth and provide detailed clues to its mechanisms of action.

## Full Story

### Issues

1. Limitations of the hard-part bio-chronology approach to reliably achieve a higher temporal resolution

2. An annual/seasonal resolution seems overly coarse, and may mask important climate signals due to excessive averaging.

Not being able to adopt a sufficiently high temporal resolution greatly constrains the power to track animal growth under climate change.

### Data

[Data policy](https://www.data.gov/privacy-policy#data_policy)

Case study tagging data [Data.gov site](https://catalog.data.gov/dataset/cooperative-shark-mark-recapture-database-mrdbs/resource/5d9ceb5f-2864-4e75-91ca-2800b402f855?inner_span=True)

[Climate indices](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/nao.shtml)

Disclaimer: Neither of the datasets are maintained by me. Please go to some official source (there are alternative data providers not listed here) to obtain your own copy of the file. Questions solely related to the dataset should be directed to the data maintainers found in the links above. However, if the links provided above are not working, let me know. The format of the current version of the data files will be different from my copy, and you may need to change some lines in the [preprocessing code](src/preprocess.R) to take that into account.

My own copy of the data files will be make available through an automatic email server by clicking [here](mailto:eidotog@gmail.com?subject=XxCLIMATE01xX&body=Do%20not%20modify%20the%20subject%20line.%20Not%20monitored.) (latency: secs). Note that I don't monitor these data requests. Once a request is fullfilled, the message will be permanantly deleted from the server. No personal information will be collected by me. There is no guarantee on the quality and timeliness of the data once it is delivered.

Datasets of this section exist in the public domain AKA CC0. They are not covered by the license of this repository. __Warning__: some data providers will try to attach unfair conditions to your data usage. Do know your rights and do not fall into those traps and choose your data provider wisely. 

### Example

The example files can be found [here](src/example). That folder contains an [R script](src/example/example.R) and a [TMB script](src/example/example.cpp). To run this example, one other [R script](src/preprocess.R), and two header files [1](src/growth.h) and [2](src/growth_imp.h) are needed.

Download the [source directory](src) first, and then make `example` you working directory. Then, you can source the [R script](src/example/example.R) to fit a von Bertalanffy curve to see if it works correctly on your machine. A verification value has been provided in the comment section at the bottom of the script. Getting exactly the same number as provided indicates success. Most likely, your R console will complain about missing libraries the first time you run it. Install them.

This code has been tested on both Windows and Linux systems. There are no system requirements to run the example code.

### Full source code

It is recommended to go through the example above first to make sure you understand the core of the code, if you are not already familiar with TMB modelling. The control statements in the full source code can be too much for the first read. But once you understand the structure of the program, the code should be fairly simple to understand.

The source files contain R scripts, [the main program](src/v3.R) and [the data preprocessor](src/preprocess.R), two header files [1](src/growth.h) and [2](src/growth_imp.h), and a [TMB script](src/v5.cpp)[^1].

The R code has been annotated and sectioned according to RStudio style. It is recommended to use RStudio to view the code.

This code has been tested on both Windows and Linux systems with 8+ GB of RAM.

__Warning__: it is recommended to run the full program on a system with 8+ GB of RAM. Some of the large models require a substantial amount of memory to execute. Your system may freeze if you are low on available memory. Your R may crash if the TMB program is given a wrong set of inputs. You have been warned. Do save your work first.

[^1]: The state of this repo is frozen such that it provides a snapshot view of the model development stage as presented in the publication. The code, however, has been continuously updated with new features added and tested. The next version of the code will appear in a different repo. Future updates of this repo are limited to minor changes, such as typo or bug fixes.