#' Title
#'
#' @param all
#' @param fine_mapping
#'
#' @return
#' @export
#'
#' @examples
get_fine_mapping_id<- function(all,fine_mapping){
CHR.all <- all$CHR
position.all <- all$position
idx.cut <- NULL
known.flag <- NULL
CHR <- fine_mapping$CHR
start <- fine_mapping$start
end <- fine_mapping$end

n.max <- nrow(all)
total <- 0
idx_cut <- rep(0,n.max)
known.flag <- rep(0,n.max)


for(i in 1:nrow(fine_mapping)){
  print(i)
  chr_temp <- CHR[i]
  start_temp <- start[i]
  end_temp <- end[i]
  idx <- which(CHR.all==chr_temp&position.all>=start_temp&
                 position.all<=end_temp)
  temp <- length(idx)
  temp.known.flag <- rep(i,length(idx))
  idx_cut[total+(1:temp)] <- idx
  known.flag[total+(1:temp)] <- temp.known.flag
  total <- temp+total
}
idx_cut <- idx_cut[1:total]
idx_cut_du <- which(duplicated(idx_cut)==F)
idx_cut <- idx_cut[idx_cut_du]
known.flag <- known.flag[1:total]
known.flag <- known.flag[idx_cut_du]
return(list(idx_cut,known.flag))

}
