#######################################
#
# Title: BioShapes
# Maintainer: L. Mada
#
# https://github.com/discoleo/BioShapes
#
# Continuation of:
# 1. Bachelor Thesis: Adrian Cotoc (2022-2023)
# Faculty of Mathematics and Informatics, UVT
#
# Coordinator:
#   Prof. Daniela Zaharie
#   Dr. med. Leonard Mada (Syonic SRL)
#   in collaboration with Syonic SRL
#   [old] GitHub: https://github.com/Adi131313/BioShapes
#
# 2. Bachelor Thesis: Darian Voda (2021-2022)



#' @example man/examples/Examples.Man.R


################
#### Arrows ####

examples.arrows()


##################
### Bio-Shapes ###

example.bioshapes()


### Liposomes

### Description of a Liposome
diagramLiposome()

### Liposome Measurements
measureLiposome()


#######################

### Enzymes: Arrows
enzymeReaction(lbl = c("A2", "B2", "Enz2", "Inhibitor 2"), new = TRUE)
enzymeReaction(y = 6, lbl = c("A1", "B1", "Enz1", "Inhibitor 1"))


#####################

### Basic Shapes

### Star Shape ###

example.star(n=5, fill = "#B0D064")


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
