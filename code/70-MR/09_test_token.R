suppressMessages(library(httr))
tok <- readLines("data/token.txt", warn = FALSE)[1]
tok <- trimws(tok)
cat("token chars:", nchar(tok), "\n")
r <- tryCatch(
  GET("https://api.opengwas.io/api/associations",
      add_headers(`X-Api-Token` = tok),
      query = list(id = "ieu-b-4869"),
      timeout(30)),
  error = function(e) e)
if (inherits(r, "error")) {
  cat("ERR:", conditionMessage(r), "\n")
} else {
  cat("HTTP", status_code(r), "\n")
  txt <- content(r, "text", encoding = "UTF-8")
  cat(substr(txt, 1, 300), "\n")
}
