context("Just testing nameKmer functionality")

test_that("Check wether nameKmer with dna alphabet for works properly",{
  names<-nameKmer(type = "dna",k=1)
  expect_equal(names,c("A","C","G","T"))

})

test_that("Check wether nameKmer with nums alphabet for works properly",{
  names<-nameKmer(type = "num",k=2,num = 3)
  expect_equal(names,c("11","12","13","21","22","23","31","32","33"))

})
