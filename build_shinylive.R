if (!requireNamespace("shinylive", quietly = TRUE)) {
  stop("Önce install.packages('shinylive') gerekli.")
}

shinylive::export(appdir = ".", destdir = "docs")

file.copy("www/ornek_madde_verisi.csv", "docs/ornek_madde_verisi.csv", overwrite = TRUE)

cat("Derleme tamam. docs/ klasörünü commit + push etmeyi unutma.\n")
