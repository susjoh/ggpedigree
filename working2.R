

library(kinship2)
library(ggplot2)
library(grid)
library(reshape)
library(plyr)

load("C:/Users/Susan Johnston/Desktop/github/ggpedigree/test_data1.Rdata")

pedigree$Sex[which(pedigree$Sex == 2)] <- "M"
pedigree$Sex[which(pedigree$Sex == 1)] <- "F"

ped <- pedigree

#ped <- read.table("../../test_78.ped", header = T)
head(ped)

ped.name.rules <- function(){
  writeLines("Pedigree object should contain the following columns:
             ID should be named ID or ANIMAL
             Mother should be MOTHER, MUM, MOM or DAM
             Father should be FATHER, DAD, POP or SIRE")
}

sex.name.rules <- function(){
  writeLines("Sex should be defined as follows:
              Females should be F or Female
              Males should be M or Male,
              Unknown should be NA, U, -9, Unk")
}


# family <- NULL


# arguments
remove.singletons = TRUE
cohort = NULL #ped$BirthYear
sex = NULL #ped$Sex
id.labels = FALSE
parent.line.key = TRUE
print.labels = TRUE
plot.unk.cohort = TRUE
single.cohort.x.shuffle = TRUE
randomise.x = TRUE
line.col.mother = "#E41A1C"
line.col.father = "#377EB8"
line.alpha = 0.3
point.size = 2
point.colour = "black"
point.alpha = 1
xlab <- ""
ylab <- ""
bg.colour = "ivory"
plot.margin = c(1.5, 1.5, 1.5, 1.5)
axis.text.colour = "darkgrey"




#~~ Format the pedigree to have ID, MOTHER, FATHER columns.

pednamevec <- c("ID", "ANIMAL", "MUM", "MOM", "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE")

names(ped)[which(toupper(names(ped)) %in% pednamevec)] <- toupper(names(ped)[which(toupper(names(ped)) %in% pednamevec)])

if(!any(c("ID", "ANIMAL")                 %in% names(ped))) stop(ped.name.rules())
if(!any(c("MUM", "MOM", "MOTHER", "DAM")  %in% names(ped))) stop(ped.name.rules())
if(!any(c("DAD", "POP", "FATHER", "SIRE") %in% names(ped))) stop(ped.name.rules())

names(ped)[which(names(ped) %in% c("ID", "ANIMAL"))]                 <- "ID"
names(ped)[which(names(ped) %in% c("MUM", "MOM", "MOTHER", "DAM"))]  <- "MOTHER"
names(ped)[which(names(ped) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"

for(i in which(names(ped) %in% c("ID", "MOTHER", "FATHER"))) ped[which(is.na(ped[,i])),i] <- 0

baseped <- ped[,c("ID", "MOTHER", "FATHER")]

#~~ Add cohort to data frame if specified

if(!is.null(cohort)) ped$graphCohort <- cohort

#~~ Add in parents that are not in the ID part as founders

baseped <- ped[,c("ID", "MOTHER", "FATHER")]

if(any(!baseped$FATHER %in% baseped$ID)) dadtab <- data.frame(ID = baseped[which(!baseped$FATHER %in% baseped$ID),"FATHER"],
                                                              MOTHER = 0, FATHER = 0)
if(any(!baseped$MOTHER %in% baseped$ID)) mumtab <- data.frame(ID = baseped[which(!baseped$MOTHER %in% baseped$ID),"MOTHER"],
                                                              MOTHER = 0, FATHER = 0)

dadtab <- unique(subset(dadtab, ID != 0))
mumtab <- unique(subset(mumtab, ID != 0))

baseped <- rbind(dadtab, mumtab, baseped)

#~~ Determine cohorts and add to the ped data.frame if not specified

if(is.null(cohort)){
  
  cohortvec <- data.frame(ID = baseped[,1],
                          graphCohort = kindepth(baseped[,"ID"],
                                                 baseped[,"FATHER"],
                                                 baseped[,"MOTHER"]))
  
}


if(!is.null(cohort)) cohortvec <- ped[,c("ID", "graphCohort")]


#~~ Remove Singletons

if(remove.singletons == TRUE){
  
  singleton.vec <- which(baseped$MOTHER == 0 & baseped$FATHER == 0 & !baseped$ID %in% c(baseped$MOTHER, baseped$FATHER))
  
  if(length(singleton.vec) > 0) baseped <- baseped[-singleton.vec,]
  
}

#~~ Melt the basepedigree and get rid of connections where value = 0

baseped2 <- melt(baseped, id.vars = c("ID"), measure.vars = c("MOTHER", "FATHER"))

baseped2 <- subset(baseped2, value != 0)

names(baseped2)[which(names(baseped2) == "value")] <- "Parent.ID"

#~~ Create a group vector for parent/offspring relationship

baseped2$Group <- 1:nrow(baseped2)
head(baseped2)

#~~ Melt to create a single line per ID with Group specified

baseped3 <- melt(baseped2, id.vars = "Group", measure.vars=c("ID", "Parent.ID"))
baseped3[1:10,]

names(baseped3)[which(names(baseped3) == "value")] <- "ID"

#~~ Recode the sex information

if(!is.null(sex)){
  sextab <- data.frame(ID = ped$ID, graphSex = sex)
  head(sextab)
  sextab$graphSex <- toupper(sextab$graphSex)
  
  sextab$graphSex[which(sextab$graphSex %in% c("M", "MALE"))]  <- 3
  sextab$graphSex[which(sextab$graphSex %in% c("F", "FEMALE"))]  <- 1
  sextab$graphSex[which(sextab$graphSex %in% c("UNKNOWN", "U", -9))]  <- 2
  sextab$graphSex[which(is.na(sextab$graphSex))] <- 2
}

#~~ Add cohort and sex information

baseped3 <- join(baseped3, cohortvec)
if(!is.null(sex)){
  baseped3 <- join(baseped3, sextab)
  baseped3$graphSex[which(is.na(baseped3$graphSex))] <- 2
}

if(is.null(sex)) baseped3$graphSex <- 1

#~~ Add parental colour information

baseped2 <- unique(subset(baseped2, select = c(Group, variable)))
names(baseped2) <- c("Group", "graphParent")
baseped3 <- join(baseped3, baseped2)


#~~ Generate X coordinates

generateXcoord <- function(size, range = c(0, 1)){
  
  if(size %% 2 != 0 & size != 1){   # Check if size is odd
    newsize <- size - 1
    interval <- diff(range)/newsize
    x <- seq(range[1], range[2], interval)
  }
  
  if(size %% 2 == 0){    # Check if size is even
    interval <- diff(range)/size
    x <- seq(range[1], range[2], interval)[-size-1] + diff(seq(range[1], range[2], interval))/2   
  }
  
  if(size == 1) x <- 0.5
  
  x
}

xcoords <- NULL

for(i in unique(baseped3$graphCohort)){
  
  # Extract the number of Unique IDs per cohort and generate X coords
  ids  <- unique(baseped3$ID[which(baseped3$graphCohort == i)])
  newx <- generateXcoord(length(ids))
  
  # Append to xcoords
  if(randomise.x == TRUE){
    xcoords <- rbind(xcoords,
                     data.frame(ID = ids,
                                xpos = sample(newx, size = length(newx), replace = F),
                                graphCohort = i))
  } else {
    xcoords <- rbind(xcoords,
                     data.frame(ID = ids,
                                xpos = newx,
                                graphCohort = i))
    
    rm(ids, newx)
  }
  
  
  if(is.na(i)){
    # Extract the number of Unique IDs per cohort and generate X coords
    ids  <- unique(baseped3$ID[which(is.na(baseped3$graphCohort))])
    newx <- generateXcoord(length(ids))
    
    
    if(randomise.x == TRUE){      
      
      # Append to xcoords
      xcoords <- rbind(xcoords,
                       data.frame(ID = ids,
                                  xpos = sample(newx, size = length(newx), replace = F),
                                  graphCohort = i))
    } else {
      
      xcoords <- rbind(xcoords,
                       data.frame(ID = ids,
                                  xpos = newx,
                                  graphCohort = i))
    }
    
    
    rm(ids, newx)
  }
  
}

# Merge with baseped3

baseped3 <- join(baseped3, xcoords)

head(baseped3)

if(single.cohort.x.shuffle == TRUE){
  
  xedit <- data.frame(Count = tapply(xcoords$graphCohort, xcoords$graphCohort, length))
  xedit$graphCohort <- row.names(xedit)
  xedit <- subset(xedit, Count == 1)
  xedit$xpos <- generateXcoord(size = nrow(xedit), range = c(0.25, 0.75))
  xedit
  
  for(i in 1:nrow(xedit)){
    baseped3[which(baseped3$graphCohort == xedit$graphCohort[i]),"xpos"] <- xedit$xpos[i]
  }
  
}

#~~ Plot

cohort.order  <- sort(unique(baseped3$graphCohort))
cohort.labels <- min(cohort.order):max(cohort.order)

if(plot.unk.cohort == TRUE){
  baseped3$graphCohort[which(is.na(baseped3$graphCohort))] <- min(cohort.order) - 1
  
  if(length(which(is.na(baseped3$graphCohort))) > 0){
    cohort.order <- c(min(cohort.order) - 1, cohort.order)
    cohort.labels <- c("Unknown", cohort.labels)
  }
}

if(print.labels == FALSE) cohort.labels <- rep("", length(cohort.order))

#~~ Recode colour if required


baseped3$graphColour <- "black"
if(parent.line.key == TRUE){
  baseped3$graphColour[which(baseped3$graphParent == "MOTHER")] <- line.col.mother
  baseped3$graphColour[which(baseped3$graphParent == "FATHER")] <- line.col.father
}  

#~~ other definitions


baseped4 <- droplevels(unique(subset(baseped3, select = c(ID, graphSex, xpos, graphCohort))))

if(id.labels == FALSE){
  ggplot() +
    geom_line(data = baseped3, aes(x = xpos, y = -graphCohort, group = Group, colour = graphColour), alpha = line.alpha) +
    geom_point(data = baseped4, aes(x = xpos, y = -graphCohort, shape = factor(graphSex)),
               size = point.size, colour = point.colour, alpha = point.alpha) +
    theme(axis.text.x      = element_blank(),
          axis.text.y      = element_text(colour = axis.text.colour),
          axis.ticks.y     = element_blank(),
          axis.ticks.x     = element_blank(),
          panel.grid       = element_blank(),
          plot.background  = element_rect(fill = bg.colour),
          panel.background = element_blank(),
          legend.position = "none") +
    theme(plot.margin = unit(plot.margin, units = "cm")) +
    scale_y_continuous(breaks = -seq(min(cohort.order), max(cohort.order), 1),
                       labels =  cohort.labels) +
    scale_colour_identity() +
    labs(x = xlab, y = ylab)
}

if(id.labels == TRUE){
  ggplot() +
    geom_line(data = baseped3, aes(x = xpos, y = -graphCohort, group = Group, colour = graphColour), alpha = line.alpha) +
    geom_text(data = baseped4, aes(x = xpos, y = -graphCohort, label = ID),
              size = point.size, colour = point.colour, alpha = point.alpha) +
    theme(axis.text.x      = element_blank(),
          axis.text.y      = element_text(colour = axis.text.colour),
          axis.ticks.y     = element_blank(),
          axis.ticks.x     = element_blank(),
          panel.grid       = element_blank(),
          plot.background  = element_rect(fill = bg.colour),
          panel.background = element_blank()) +
    theme(plot.margin = unit(plot.margin, units = "cm")) +
    scale_y_continuous(breaks = -seq(min(cohort.order), max(cohort.order), 1),
                       labels =  cohort.labels) +
    scale_colour_identity() +
    labs(x = xlab, y = ylab)
}
