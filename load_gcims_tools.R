# Load all custom GC-IMS functions from the R/ folder
library(here)
cat(">>> loading GCIMS tools...\n")
files <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
print(files)
for (f in files) {
  cat("→ Sourcing:", f, "\n")
  tryCatch(
    source(f),
    error = function(e) {
      cat("⚠️ Error sourcing", f, ":", e$message, "\n")
    }
  )
}
