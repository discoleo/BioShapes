###################
#
# Bachelor Thesis
#
# Title: BioShapes
#
# Candidate: Adrian Cotoc
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#
# in collaboration with Syonic SRL
# continous the work of Darian Voda
#
# GitHub: https://github.com/Adi131313/BioShapes


#' @example man/examples/Examples.Man.R


################
#### Arrows ####

examples.arrows()


####################
#### Bio-Shapes ####

examples.bioshapes()


### Liposomes

### Description of a Liposome

diagramLiposome()

### Liposome Measurements

measureLiposome()


#######################

### Enzymes: Arrows
plot.base()
enzymeReaction(lbl = c("A2", "B2", "Enz2", "Inhibitor 2"))
enzymeReaction(y = 6, lbl = c("A1", "B1", "Enz1", "Inhibitor 1"))


#####################

### Basic Shapes

### Star Shape ###

example.star(n=5, fill = "yellow")


### Arcs

### Simple Arcs
example.arcs()

### Apple-Melon Inversion:
# (weird things can happen at 1 MT)
example.ArcsByDist(lwd=2)


### Lens Examples ###

example.lens()



### Various Curves ###

example.curves()


### Simple Braces ###

example.braces()


#####################

### Objects

### DNA Structure ###

example.dna()


#####################

### Neuron

### Creation of the Neuron ###

example.neuronDesign()

### Neuron: Basic ###

example.neuron()

### Description of the Neuron ###

description.neuron()

### Multiple neurons ###
# TODO: move neurons around and make connections between
example.neurons()


### Muscle Tissue ###

example.muscle()


### Blood Cells ###

example.bloodCells()


### Duct Examples ###

example.ducts()

### Complex Duct ###

example.complexDuct(n = 8)

### Viruses ###

example.virus()

#####################
