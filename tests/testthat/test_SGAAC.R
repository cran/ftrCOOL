context("Just testing Seprated Grouped Amino Acid composition functionality")

test_that("Check whether SGAAC works properly",{

  comp_part2<-as.vector(SGAAC("ACYMELACYL",numNterm = 2,numCterm = 2,normalized = FALSE,Grp = "aromatic"))
  #1521411521
  #ACYMELACYL
  excomp_p1<-rep(0,5)
  excomp_p2<-rep(0,5)
  excomp_p3<-rep(0,5)
  namesP<-c("1","2","3","4","5")
  names(excomp_p1)<-namesP
  names(excomp_p2)<-namesP
  names(excomp_p3)<-namesP

  excomp_p1[c("1","5")]=c(1,1)
  excomp_p3[c("2","1")]=c(1,1)
  excomp_p2[c("2","1","4","5")]=c(1,3,1,1)

  expected<-c(excomp_p1,excomp_p2,excomp_p3)
  names(expected)<-NULL

  expect_equal(comp_part2,expected)

  expect_error(SGAAC("ACYMELACYL",numNterm = 10,numCterm = 2,normalized = FALSE,Grp="aromatic"))

  expect_warning(SGAAC(c("ACYMELACYL","ACHTCILNMSAHKGHCV"),numNterm = 10,numCterm = 2,normalized = FALSE,Grp="aromatic"))


})
