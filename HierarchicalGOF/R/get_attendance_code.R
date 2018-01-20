#' @title Open code and run the attendance analysis
#' @author Devin S. Johnson
#' @export

run_attendance_analysis = function(show=TRUE, run=FALSE){
  
  dir = system.file("CSL_attendance", package="HierarchicalGOF")
  
  if(show) file.edit(file.path(dir, "find_post_mode_hsmm.R"))
  if(run) source(file.path(dir, "find_post_mode_hsmm.R"), local=TRUE)
  
  if(show) file.edit(file.path(dir, "sample_post_hsmm.R"))
  if(run) source(file.path(dir, "sample_post_hsmm.R"), local=TRUE)
  
  if(show) file.edit(file.path(dir, "attendance_plots.R"), local=TRUE)
  if(run) source(file.path(dir, "attendance_plots.R"))
  if(run) system(paste0('open "', file.path(dir, "attendance_fit.pdf"), '"'))
  
}
