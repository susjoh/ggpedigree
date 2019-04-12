#' ggpedigree: Plot a pedigree
#' @param pedigree A data frame of id, dam and sire
#' @param cohort Default NULL. An optional vector assigning a cohort to each id.
#'   If NULL then uses kindepth function in kinship2 to assign IDs to cohorts.
#' @param sex Default NULL. An optional vector assigning a sex to each id
#' @param family Default NULL. An optional vector assigning a family to each id
#'   (NOT WORKING YET)
#' @param id.labels logical. Default FALSE. Print the IDs on the pedigree.
#' @param remove.singletons logical. Default TRUE. Remove IDs with no offspring
#'   or parents assigned.
#' @param colour.parent.lines logical. Default TRUE. Colours lines based on
#'   maternal and paternal links.
#' @param print.cohort.labels logical. Default TRUE. Prints cohort ideas on left
#'   hand side of plot.
#' @param plot.unk.cohort logical. Default TRUE. Show IDs of unknown cohort.
#' @param single.cohort.x.shuffle logical. Default TRUE. Randomly assign left to
#'   right position of IDs.
#' @param randomise.x logical. Default TRUE.
#' @param return.plot.tables logical. Default FALSE. Returns an object with the
#'   line and point data used for the plot.
#' @param suppress.plot logical. Default FALSE. Print plot or not.
#' @param line.col.mother line colour for mother
#' @param line.col.father line colour for father
#' @param line.alpha line alpha for plot
#' @param point.size  point size
#' @param point.colour point colour
#' @param point.alpha point alpha
#' @param gg.theme extra theme information for ggplot
#' @param xlab x axis label
#' @param ylab y axis label
#' @param bg.colour background colour
#' @param plot.margin plot margin dimensions
#' @param axis.text.colour axis text colour
#' @import reshape2
#' @export


ggpedigree <- function(pedigree, 
                       cohort                  = NULL, 
                       sex                     = NULL,
                       family                  = NULL,
                       id.labels               = FALSE,
                       remove.singletons       = TRUE,
                       colour.parent.lines     = TRUE,
                       print.cohort.labels     = TRUE,
                       plot.unk.cohort         = TRUE,
                       single.cohort.x.shuffle = TRUE,
                       randomise.x             = TRUE,
                       return.plot.tables      = FALSE,
                       suppress.plot           = FALSE,
                       line.col.mother         = "#E41A1C",
                       line.col.father         = "#377EB8",
                       line.alpha              = 0.3,
                       point.size              = 2,
                       point.colour            = "black",
                       point.alpha             = 1,
                       gg.theme                = NULL,
                       xlab                    = "",
                       ylab                    = "",
                       bg.colour               = "ivory",
                       plot.margin             = c(1.5, 1.5, 1.5, 1.5),
                       axis.text.colour        = "darkgrey"){
   
  ped <- pedigree
  
  #~~ Format the pedigree to have ID, MOTHER, FATHER columns and recode NA to 0.
  
  pednamevec <- c("ID", "ANIMAL", "MUM", "MOM", "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE")
  
  names(ped)[which(toupper(names(ped)) %in% pednamevec)] <- toupper(names(ped)[which(toupper(names(ped)) %in% pednamevec)])
  
  
  if(!any(c("ID", "ANIMAL")                 %in% names(ped))) stop(ped.name.rules())
  if(!any(c("MUM", "MOM", "MOTHER", "DAM")  %in% names(ped))) stop(ped.name.rules())
  if(!any(c("DAD", "POP", "FATHER", "SIRE") %in% names(ped))) stop(ped.name.rules())
  
  names(ped)[which(names(ped) %in% c("ID", "ANIMAL"))]                 <- "ID"
  names(ped)[which(names(ped) %in% c("MUM", "MOM", "MOTHER", "DAM"))]  <- "MOTHER"
  names(ped)[which(names(ped) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"
  
  for(i in which(names(ped) %in% c("ID", "MOTHER", "FATHER"))) ped[,i] <- as.character(ped[,i])
  for(i in which(names(ped) %in% c("ID", "MOTHER", "FATHER"))) ped[which(is.na(ped[,i])),i] <- 0
  
  #~~ Check that ID has not been duplicated
  
  if( is.null(family)) if(any(as.numeric(names(table(table(ped$ID)))) > 1)) stop("Duplicated values in ID column")
  if(!is.null(family)) if(any(as.numeric(names(table(table(paste(ped$Family, ped$ID))))) > 1)) stop("Duplicated values in ID column")
  
  #~~ Create a baseped object
  
  baseped <- ped[,c("ID", "MOTHER", "FATHER")]
  baseped$MOTHER[is.na(baseped$MOTHER)] <- 0
  baseped$FATHER[is.na(baseped$FATHER)] <- 0
  
  #~~ Add in parents that are not in the ID part as founders
  
  if(any(!baseped$FATHER %in% baseped$ID)) dadtab <- data.frame(ID = baseped[which(!baseped$FATHER %in% baseped$ID),"FATHER"],
                                                                MOTHER = 0, FATHER = 0)
  if(any(!baseped$MOTHER %in% baseped$ID)) mumtab <- data.frame(ID = baseped[which(!baseped$MOTHER %in% baseped$ID),"MOTHER"],
                                                                MOTHER = 0, FATHER = 0)
  
  dadtab <- unique(subset(dadtab, ID != 0))
  mumtab <- unique(subset(mumtab, ID != 0))
  
  baseped <- rbind(dadtab, mumtab, baseped)
  rm(dadtab, mumtab)
  
  #~~ Determine cohorts and add to the ped data.frame if not specified
  
  if(is.null(cohort)){
    
    cohortvec <- data.frame(ID = baseped[,1],
                            graphCohort = kindepth(baseped[,"ID"],
                                                   baseped[,"FATHER"],
                                                   baseped[,"MOTHER"]))
    
  }
  
  
  if(!is.null(cohort)) cohortvec <- data.frame(ID = ped[,"ID"], graphCohort = cohort)
  
  #~~ Recode the sex information
  
  if(!is.null(sex)){
    sextab <- data.frame(ID = ped$ID, graphSex = sex)
    head(sextab)
    sextab$graphSex <- toupper(sextab$graphSex)
    
    if(any(!sextab$graphSex %in% c("M", "F", "MALE", "FEMALE", "UNKNOWN", "U", "-9", NA))) stop(sex.name.rules())
    
    sextab$graphSex[which(sextab$graphSex %in% c("M", "MALE"))]  <- 3
    sextab$graphSex[which(sextab$graphSex %in% c("F", "FEMALE"))]  <- 1
    sextab$graphSex[which(sextab$graphSex %in% c("UNKNOWN", "U", -9))]  <- 2
    sextab$graphSex[which(is.na(sextab$graphSex))] <- 2
  }
  
  #~~ Remove Singletons
  
  if(remove.singletons == TRUE){
    
    singleton.vec <- which(baseped$MOTHER == 0 & baseped$FATHER == 0 & !baseped$ID %in% c(baseped$MOTHER, baseped$FATHER))
    
    if(length(singleton.vec) > 0) baseped <- baseped[-singleton.vec,]
    
  }
  
  #~~ Melt baseped and get rid of connections where value = 0 (means parental connection is unknown)
  
  baseped2 <- melt(baseped, id.vars = c("ID"), measure.vars = c("MOTHER", "FATHER"))
  
  baseped2 <- subset(baseped2, value != 0)
  
  names(baseped2)[which(names(baseped2) == "value")] <- "Parent.ID"
  
  #~~ Create a group vector for parent/offspring relationship
  
  baseped2$Group <- 1:nrow(baseped2)
  
  #~~ Melt to create a single line per ID with Group specified
  
  baseped2$ID <- as.character(baseped2$ID)
  baseped3 <- melt(baseped2, id.vars = "Group", measure.vars=c("ID", "Parent.ID"))
  baseped3[1:10,]
  
  names(baseped3)[which(names(baseped3) == "value")] <- "ID"
  
  
  
  #~~ Add cohort and sex information
  
  suppressMessages(baseped3 <- join(baseped3, cohortvec))
  if(!is.null(sex)){
    suppressMessages(baseped3 <- join(baseped3, sextab))
    baseped3$graphSex[which(is.na(baseped3$graphSex))] <- 2
  }
  
  if(is.null(sex)) baseped3$graphSex <- 1
  
  #~~ Add parental colour information
  
  baseped2 <- unique(subset(baseped2, select = c(Group, variable)))
  names(baseped2) <- c("Group", "graphParent")
  suppressMessages(baseped3 <- join(baseped3, baseped2))
  
  
  baseped3$graphColour <- "black"
  if(colour.parent.lines == TRUE){
    baseped3$graphColour[which(baseped3$graphParent == "MOTHER")] <- line.col.mother
    baseped3$graphColour[which(baseped3$graphParent == "FATHER")] <- line.col.father
  }  
  
  
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
    
    if(!is.na(i)){
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
  
  #~~ Merge with baseped3
  
  suppressMessages(baseped3 <- join(baseped3, xcoords))
  
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
  
  #~~ Collate information for axis labels and plotting unknown cohort individuals
  
  cohort.order  <- sort(unique(baseped3$graphCohort))
  cohort.labels <- min(cohort.order):max(cohort.order)
  
  if(plot.unk.cohort == TRUE){
    baseped3$graphCohort[which(is.na(baseped3$graphCohort))] <- min(cohort.order) - 1
    
    if(length(which(baseped3$graphCohort == min(cohort.order) - 1)) > 0){
      cohort.order <- c(min(cohort.order) - 1, cohort.order)
      cohort.labels <- c("Unknown", cohort.labels)
    }
  }
  
  if(print.cohort.labels == FALSE) cohort.labels <- rep("", length(cohort.order))
  
  
  #~~ baseped3 is used for the parental links. Create baseped4 which is a unique value for 
  #   each individual to avoid overplotting.
  
  baseped4 <- droplevels(unique(subset(baseped3, select = c(ID, graphSex, xpos, graphCohort))))
  
  
  #~~ Create a return object if return.plot.tables == TRUE
  
  if(return.plot.tables == TRUE) return(list(LinePlotFrame = baseped3,
                                            PointPlotFrame = baseped4))
  
  #~~ Plot the pedigrees
  
  if(is.null(gg.theme)) {
    gg.theme <- theme(axis.text.x      = element_blank(),
        axis.text.y      = element_text(colour = axis.text.colour),
        axis.ticks.y     = element_blank(),
        axis.ticks.x     = element_blank(),
        panel.grid       = element_blank(),
        plot.background  = element_rect(fill = bg.colour),
        panel.background = element_blank(),
        legend.position = "none",
        plot.margin = unit(plot.margin, units = "cm"))
  }
  
  
  if(suppress.plot == FALSE){
    
    p <- ggplot() +
      geom_line(data = baseped3, aes(x = xpos, y = -graphCohort, group = Group, colour = graphColour), alpha = line.alpha) +
      scale_y_continuous(breaks = -seq(min(cohort.order), max(cohort.order), 1),
                         labels =  cohort.labels) +
      scale_colour_identity() +
      labs(x = xlab, y = ylab) +
      gg.theme
    
    
    
    if(id.labels == FALSE){
      print(p +
        geom_point(data = baseped4, aes(x = xpos, y = -graphCohort, shape = factor(graphSex)),
                   size = point.size, colour = point.colour, alpha = point.alpha))
        
    }
    
    if(id.labels == TRUE){
      print(p +
        geom_text(data = baseped4, aes(x = xpos, y = -graphCohort, label = ID),
                  size = point.size, colour = point.colour, alpha = point.alpha))
    }
  }
}



