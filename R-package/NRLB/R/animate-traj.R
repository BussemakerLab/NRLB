#' Animate solution trajectory with PSAMs
#' 
#' @param trajPath TO_BE_ADDED
#' @param outPath TO_BE_ADDED
#' @param k TO_BE_ADDED
#' @return TO_BE_ADDED
#' 
#' @examples
#'
#' @export
#' 
animate.traj = function(trajPath, outPath, k) {
  traj = read.table(trajPath, header=FALSE)
  dir.create(path="~/tmp_animate_dir")
  for (i in 1:nrow(traj)) {
    logo(betas=traj[i,1:(k*4)], save.path = paste0("~/tmp_animate_dir/",formatC(i, width=5, format="d", flag="0"),".png"), display = FALSE, isPDF = FALSE, title=i)
  }
  system(paste0("convert -delay 80 ~/tmp_animate_dir/*.png ", outPath))
  system("rm -rf ~/tmp_animate_dir")  
}
