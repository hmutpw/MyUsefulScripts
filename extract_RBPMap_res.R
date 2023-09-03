get_binding_sites <- function(filePath){
  lines <- readLines(filePath)
  
  result <- list()
  current_group <- NULL
  
  group_start_lines <- which(startsWith(lines, "==="))
  group_end_lines <- which(startsWith(lines, "***"))[-1]
  group_id <- lines[group_start_lines-1]
  group_id <- sub("^\\s*=+\\s*", "", group_id)
  group_id <- sub("^\\s+", "", group_id)
  content_lines <- group_end_lines - group_start_lines
  group_id_long <- rep(group_id,content_lines)
  
  
  
  out_df <- purrr::map2_df(.x = group_start_lines, .y = group_end_lines, .f = function(st, ed, lines){
    # get group id
    group_id <- lines[st-1]
    group_id <- sub("^\\s*=+\\s*", "", group_id)
    group_id <- sub("^\\s+", "", group_id)
    current_lines <- lines[(st+1):(ed-1)]
    current_lines <- current_lines[current_lines !=""]
    protein_name_line <- which(startsWith(current_lines, "Protein:"))
    protein_names <- sub("^Protein:\\s+(.*?)\\(Hs\\/Mm\\)", "\\1", current_lines[protein_name_line])
    protein_names_long <- rep(protein_names, diff(c(protein_name_line,length(current_lines)+1))-1)
    current_line_list <- split(current_lines[-protein_name_line],protein_names_long)
    current_line_list <- lapply(current_line_list, function(x){
      df <- stringr::str_split(string = x, pattern = "\t",simplify = T)
      df <- t(apply(df,1,trimws))
      colnames(df) <- df[1,]
      as.data.frame(df[-1,,drop=FALSE])
    })
    current_line_df <- data.table::rbindlist(l = current_line_list,idcol = "Protein")
    current_line_df$event_id <- rep(group_id,nrow(current_line_df))
    current_line_df
  },lines = lines)
}
