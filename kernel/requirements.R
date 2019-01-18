packs_to_install = c(
	'evolqg', 
	'devtools',
	'phytools',
	'geomorph',
	'shapes',
	'RPostgreSQL'
)

install.packages(packs_to_install, dependencies=TRUE, repos='http://cran.rstudio.com/')
