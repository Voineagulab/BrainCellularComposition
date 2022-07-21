
## First: check if the root directory is defined
if (!(base::exists(root.dir))) {
  stop("Please define the object \"root.dir\" to be root directory of all these replication analyses. It should contain the directory structure as described on our GitHub")
}

## Plotting
  invis <- element_blank() # a shorthand for ggplot2
  