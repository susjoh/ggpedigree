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