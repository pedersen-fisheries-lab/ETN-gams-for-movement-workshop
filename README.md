# ETN Introduction to GAMs workshop

A (very!) short course on how to fit, plot, and evaluate GAMs

The course website URL  available at [pedersen-fisheries-lab.github.io/ETN-gams-for-movement-workshop/](https://pedersen-fisheries-lab.github.io/ETN-gams-for-movement-workshop/)

## Setup

  - You will need to install **R** and I recommend using **RStudio**. The
    latest version of R can be downloaded
    [here](https://cran.r-project.org/mirrors.html). RStudio is an application
    (an integrated development environment or IDE) that facilitates the use of R
    and offers a number of nice additional features. It can be downloaded
    [here](https://www.rstudio.com/products/rstudio/download/). You will need
    the free Desktop version for your computer.

  - Download the course materials as a ZIP file
    [here](https://github.com/pedersen-fisheries-lab/ETN-gams-for-movement-workshop/archive/main.zip).
    Alternatively, if you have the [**usethis**](), R package, running the
    following command will download the course materials and open them:

    ``` {.r}
    usethis::use_course('pedersen-fisheries-lab/ETN-gams-for-movement-workshop')
    ```

  - Install the R packages required for this course by running the following
    line of code your R console:

    ``` {.r}
    install.packages(c("dplyr", "ggplot2", "sf", "mgcv", "tidyr", "gratia")
    ```
    
## Part 1: What is a GAM? Basis functions and smoothers


[Introduction slides](slides/intro.html)

[Part 1 slides](slides/part1_Intro_to_GAMs.html)

[Part 1 script](scripts/part1_intro_to_gams.R)



## Part 2

[Part 2 slides](slides/Part2_GAMs_for_movement_data.html)

[Part 2 script](scripts/part2_families_and_interactions.R)

## Extra resources:

If you are looking for more GAM teaching aids, here's a few resources that might be helpful:


### Longer workshops:

* [Older 1-day GAM workshop for the Ecological Society of America](http://eric-pedersen.github.io/mgcv-esa-workshop/)
* [Recent 1-day GAM workshop](https://pedersen-fisheries-lab.github.io/one-day-gam-workshop/)
* [Recent 3-day workshop for DFO Canada](https://github.com/pedersen-fisheries-lab/DFO-3day-gam-workshop)


### Video tutorials:

* Gavin Simpson: 
  - [Youtube workshop 1](https://www.youtube.com/watch?v=sgw4cu8hrZM)
  - [Youtube workshop 2](https://www.youtube.com/watch?v=Ukfvd8akfco)
    
### Interactive tutorials:

* Noam Ross:
  - [GAMS in R Course](https://noamross.github.io/gams-in-r-course/)
    
### Written resources:

* [David Lawrence Miller's blog](https://converged.yt/)
* [Gavin Simpson's blog](https://fromthebottomoftheheap.net/)
* Simon Wood (The creator of mgcv): papers and books
  - "Generalized Additive Models: An Introduction with R 2nd ed"
  - 2020: TEST: "Inference and computation with generalized additive models and their extensions"
