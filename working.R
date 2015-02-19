

library(kinship2)
library(ggplot2)
library(grid)
library(reshape)
library(plyr)

load("C:/Users/Susan Johnston/Desktop/github/ggpedigree/test_data1.Rdata")

# ped <- pedigree
ped <- pedigree

simple.ped.name.rules <- function(){
  writeLines("Pedigree object should contain the following columns:
               ID should be named ID or ANIMAL
               Mother should be MOTHER, MUM, MOM or DAM
               Father should be FATHER, DAD, POP or SIRE")
}

# arguments
remove.singletons <- TRUE
cohort <- NULL
sex <- ped$Sex
family <- ped$Family
print.labels <- TRUE

#~~ Format the pedigree

pednamevec <- c("ID", "ANIMAL", "MUM", "MOM", "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE")

names(ped)[which(toupper(names(ped)) %in% pednamevec)] <- toupper(names(ped)[which(toupper(names(ped)) %in% pednamevec)])

if(!any(c("ID", "ANIMAL") %in% names(ped)))   stop(simple.ped.name.rules())
if(!any(c("MUM", "MOM", "MOTHER", "DAM") %in% names(ped))) stop(simple.ped.name.rules())
if(!any(c("DAD", "POP", "FATHER", "SIRE") %in% names(ped))) stop(simple.ped.name.rules())

names(ped)[which(names(ped) %in% c("ID", "ANIMAL"))]                 <- "ID"
names(ped)[which(names(ped) %in% c("MUM", "MOM", "MOTHER", "DAM"))]  <- "MOTHER"
names(ped)[which(names(ped) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"

for(i in which(names(ped) %in% c("ID", "MOTHER", "FATHER"))) ped[which(is.na(ped[,i])),i] <- 0

#~~ Determine cohorts and add to the ped data.frame

if(is.null(cohort)){
  
  cohortvec <- data.frame(ID = ped[,1],
                          graphCohort = kindepth(ped[,"ID"], ped[,"FATHER"], ped[,"MOTHER"]))
  
  ped <- join(ped, cohortvec)
  
}

if(!is.null(cohort)) ped$graphCohort <- cohort




#~~ Add sex to data frame if specified

if(!is.null(sex)) ped$graphSex <- sex

#~~ Remove Singletons

if(remove.singletons == TRUE){
  singleton.vec <- which(ped$MOTHER == 0 & ped$FATHER == 0 & !ped$ID %in% c(ped$MOTHER, ped$FATHER))
  if(length(singleton.vec) > 0) ped <- ped[-singleton.vec,]
}

#~~ Melt the pedigree

ped2 <- melt(ped, id.vars = c("ID"), measure.vars = c("MOTHER", "FATHER"))
ped2 <- subset(ped2, !is.na(value))

head(ped2)

names(ped2)[which(names(ped2) == "value")] <- "Parent.ID"
ped2$Group <- 1:nrow(ped2)
head(ped2)

ped3 <- melt(ped2, id.vars = "Group", measure.vars=c("ID", "Parent.ID"))
ped3[1:10,]

names(ped3)[3] <- "ID"

#~~ Add cohort and sex information

ped3 <- join(ped3, ped[,c("ID", "graphCohort")])

if(!is.null(sex)) ped3 <- join(ped3, ped[,c("ID", "graphSex")])

head(ped3)

#~~ Remove groups with zeros

if(0 %in% ped3$ID){
  groupRemove <- ped3[which(ped3$ID == 0), "Group"]
  ped3 <- subset(ped3, !Group %in% groupRemove)
}

#~~ Generate X coordinates

generateXcoord <- function(size){
  
  if(size %% 2 != 0 & size != 1){   # Check if size is odd
    newsize <- size - 1
    interval <- 1/newsize
    x <- seq(0, 1, interval)
  }
  
  if(size %% 2 == 0){    # Check if size is even
    interval <- 1/size
    x <- seq(0, 1, interval)[-size-1] + diff(seq(0, 1, interval))/2   
  }
  
  if(size == 1) x <- 0.5
  
  x
}


xcoords <- NULL

for(i in unique(ped3$graphCohort)){
  
  # Extract the number of Unique IDs per year and generate X coords
  ids  <- unique(ped3$ID[which(ped3$graphCohort == i)])
  newx <- generateXcoord(length(ids)) # generate X coordinates
  
  # Append to xcoords
  xcoords <- rbind(xcoords,
                   data.frame(ID = ids,
                              x = sample(newx, size = length(newx), replace = F)))
  
  rm(ids, newx)
}

# Merge with ped3

ped3 <- join(ped3, xcoords)

#~~ Plot

cohort.order  <- sort(unique(ped3$graphCohort))
cohort.labels <- cohort.order

if(print.labels == FALSE) cohort.labels <- rep("", length(cohort.order))

line.alpha <- 0.1
point.size <- 4
xlab <- ""
ylab <- ""
legend.pos <- ifelse(is.null(sex), "none", "bottom")
bg.colour <- "ivory"

ggplot(ped3, aes(x, -graphCohort)) +
  geom_line(aes(group = Group), alpha = line.alpha) +
  geom_point(size = point.size) +
  theme(axis.text.x      = element_blank(),
        axis.text.y      = element_text(colour = "darkgrey"),
        axis.ticks.y     = element_blank(),
        axis.ticks.x     = element_blank(),
        panel.grid       = element_blank(),
        plot.background  = element_rect(fill = bg.colour),
        panel.background = element_blank()) +
  theme(plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), units = "cm")) +
  scale_y_continuous(breaks = -seq(min(cohort.order), max(cohort.order), 1),
                     labels =  cohort.labels) +
  labs(x = xlab, y = ylab)









##########################################


} else {
  
  newped <- NULL
  
  for(i in unique(family)){
    
    tempped <- subset(ped, family == i)
    
    cohortvec <- data.frame(ID = tempped[,1],
                            graphCohort = kindepth(tempped[,"ID"],
                                                   tempped[,"FATHER"],
                                                   tempped[,"MOTHER"]))
    
    tempped <- join(tempped, cohortvec)
    
    newped <- rbind(newped, tempped)
    
    rm(tempped, cohortvec)
  }
}

newped
