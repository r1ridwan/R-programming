library(GenomicRanges)
library(plyranges)

## Create a GRanges object
gr <- GRanges(seqnames = "chr1", 
              strand = c("+", "-"),
              ranges = IRanges(start = c(102012,520211),
                               end=c(120303, 526211)), 
              gene_id = c(1001,2151), 
              score = c(10, 25))
gr


## Create a new IRanges object
ir_new <- IRanges(start = c(1:3),end=c(5,2,8))
ir_new

## Modify (or mutate) a genomic interval by altering the width
## Use the mutate() verb with anchor_*() adverbs
mutate(ir_new, width = 10) # width will be 10
mutate(anchor_start(ir_new), width = 10) # end value will change to make width 10
mutate(anchor_end(ir_new), width = 10) # start value will change to make width 10

## Anchoring fixes either the start, end, center coordinates 
##    of the GRange object. 
## Can be used in combination with either mutate or stretch


## Stretch a genomic interval by an intger amount 
ir_new
stretch(anchor_center(ir_new), extend = 10)

## Shift (left and right)
ir_new
shift_left(ir_new, shift = 10)
shift_right(ir_new, shift = 10)


## Create a GRanges object
gr <- 
  data.frame(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE),
             strand = sample(c("+", "-"), 7, replace = TRUE),
             score = runif(7),
             start = 1:7,
             width = 10) %>%
  as_granges()
gr

## Summarizing by group
gr %>%
  group_by(strand) %>% 
  summarize(mean_score = mean(score))

## Filtering by group
gr %>%
  group_by(strand) %>% 
  filter(score == max(score))
