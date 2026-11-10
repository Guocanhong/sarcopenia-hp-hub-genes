suppressMessages(library(httr))
tok <- trimws(readLines("data/token.txt", warn = FALSE)[1])
cat("token head:", substr(tok, 1, 25), "...\n")
for (ep in c("https://api.opengwas.io/api/status",
             "https://api.opengwas.io/api/user",
             "https://api.opengwas.io/api/gwasinfo/ieu-b-4869")) {
  r <- tryCatch(GET(ep, add_headers(`X-Api-Token` = tok), timeout(30)),
                error = function(e) e)
  if (inherits(r, "error")) { cat(ep, "ERR:", conditionMessage(r), "\n"); next }
  txt <- content(r, "text", encoding = "UTF-8")
  cat(ep, "-> HTTP", status_code(r), ":", substr(txt, 1, 150), "\n")
}
