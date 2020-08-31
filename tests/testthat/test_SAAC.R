context("Just testing Seprated Amino Acid composition functionality")

test_that("Check whether SAAC works properly",{

  comp_part2<-as.vector(SAAC("ACYMELACYL",numNterm = 2,numCterm = 2,normalized = FALSE))

  excomp_p1<-rep(0,20)
  excomp_p2<-rep(0,20)
  excomp_p3<-rep(0,20)
  namesP<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  names(excomp_p1)<-namesP
  names(excomp_p2)<-namesP
  names(excomp_p3)<-namesP

  excomp_p1[c("A","C")]=c(1,1)
  excomp_p3[c("Y","L")]=c(1,1)
  excomp_p2[c("Y","M","E","L","A","C")]=c(1,1,1,1,1,1)

  expected<-c(excomp_p1,excomp_p2,excomp_p3)
  names(expected)<-NULL

  expect_equal(comp_part2,expected)

  expect_error(SAAC("ACYMELACYL",numNterm = 10,numCterm = 2,normalized = FALSE))

  expect_warning(SAAC(c("ACYMELACYL","ACHTCILNMSAHKGHCV"),numNterm = 10,numCterm = 2,normalized = FALSE))


})
