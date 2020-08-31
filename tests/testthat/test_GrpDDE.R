context("Just testing GrpDDE functionality")

test_that("Check whether GrpDDE works properly",{

  dde<-GrpDDE(seqs="AYAAA",Grp = "aromatic")
  #12111

  #AY=4*2/61 YA=2*4/61 AA=4*4/61 AA=4*4/61
  #("A"=4,"C"=2,"D"=2,"E"=2,"F"=2,"G"=4,"H"=2,"I"=3,"K"=2,"L"=6,"M"=1,"N"=2,"P"=4,"Q"=2,"R"=6,"S"=6,"T"=4,"V"=4,"W"=1,"Y"=2)


  #Grp1=c("G","A","V","L","M","I") 4+4+4+6+1+3 = 22
  #,Grp2=c("F","Y","W") 2+2+1 =5

  #12= 22*5      21=  5*22       11=22*22       11= 22*22

  tm_11=(22*22)/(61*61)
  dc_11=2/4
  tv_11=(tm_11)*(1-tm_11)/(4)
  dde_11=(dc_11-tm_11)/sqrt(tv_11)

  expect_equal(dde[1,"11"],dde_11)



})
