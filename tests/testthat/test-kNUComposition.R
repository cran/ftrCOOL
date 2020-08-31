context("Just testing kNUComposition functionality")

test_that("Check whether NUComposition works properly k=1",{

  knucompos<-as.vector(kNUComposition(seqs="AGCCCCGT",rng = 1,normalized = FALSE))
  expected<-vector(mode = "numeric",length = 4)
  names(expected)<-nameKmer(k=1,type = "dna")
  expected[c("A","C","G","T")]=c(1,4,2,1)
  names(expected)<-NULL
  expect_equal(knucompos,expected)
})

test_that("Check whether DiNucleotide Composition works properly",{
  knucompos<-as.vector(kNUComposition(seqs="AGCCCCGT",rng = 2,normalized = FALSE))
  expected<-vector(mode = "numeric",length = 16)
  names(expected)<-nameKmer(k=2,type = "dna")
  expected[c("AG","GC","CC","CG","GT")]=c(1,1,3,1,1)
  names(expected)<-NULL
  expect_equal(knucompos,expected)
})
