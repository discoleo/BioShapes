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



example.arrows()


####################
#### Bio-Shapes ####

examples.bioshapes()


#### Description of a Liposome ####

diagramLiposome()

#### Liposome Measurements ####

measureLiposome()


### Arrows with Enzymes ####
plot.base()
enzymeReaction(lbl = c("A2", "B2", "Enz2", "Inhibitor 2"))
enzymeReaction(y = 6, lbl = c("A1", "B1", "Enz1", "Inhibitor 1"))

### DNA Structure ###

example.dna()

### Blood Cells ###

example.bloodCells()

### Creation of the neuron ###

example.neuronDesign()

### Neuron ###

example.neuron()

### Description of the Neuron ###

description.neuron()

### Multiple neurons ###
# TODO: move neurons around and make connections between
example.neurons()

### Muscle Tissue ###

example.muscle()

### Star Shape ###

example.star(n=5, fill = "yellow")

### Simple Braces ###

example.braces()

### Various Curves ###

example.curves()

### Arcs Examples ###
# TODO: create the 4th lens
example.arcs()

### Lens Examples ###

example.lens()

### Duct Examples ###

example.ducts()

### Complex Duct ###

example.complexDuct(n = 8)

### Viruses ###

example.virus()

#####################
