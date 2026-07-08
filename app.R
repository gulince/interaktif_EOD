
library(shiny)
library(ggplot2)
library(bslib)

ham_alpha_hesapla <- function(df_sayisal) {
  df_sayisal <- as.matrix(df_sayisal)
  k <- ncol(df_sayisal)
  kovaryans <- cov(df_sayisal, use = "pairwise.complete.obs")
  toplam_varyans <- sum(kovaryans)
  madde_varyans_toplami <- sum(diag(kovaryans))
  (k / (k - 1)) * (1 - madde_varyans_toplami / toplam_varyans)
}

dagilim_uret <- function(n, carpiklik, basiklik_duzey, hedef_ortalama, hedef_sd) {
  if (basiklik_duzey > 0) {
    taban <- rt(n, df = 60 - basiklik_duzey * 5.4)
  } else {
    # Beta(p,p) (p = 0.8076923) tam olarak -1.3 fazla basiklik verir;
    # normal (basiklik=0) ile karistirilarak duzey -10..0 araligi kapsanir.
    agirlik <- -basiklik_duzey / 10
    p <- 0.8076923
    sd_beta <- 1 / (2 * sqrt(2 * p + 1))
    bilesen <- (rbeta(n, p, p) - 0.5) / sd_beta
    taban <- (1 - agirlik) * rnorm(n) + agirlik * bilesen
  }
  olcekli <- ifelse(taban < 0, taban * (1 - carpiklik), taban * (1 + carpiklik))
  standart <- (olcekli - mean(olcekli)) / sd(olcekli)
  standart * hedef_sd + hedef_ortalama
}

mod_hesapla <- function(x) {
  yogunluk <- density(x)
  yogunluk$x[which.max(yogunluk$y)]
}

skewness_hesapla <- function(x) {
  n <- length(x)
  ortalama <- mean(x)
  s <- sqrt(sum((x - ortalama)^2) / n)
  (sum((x - ortalama)^3) / n) / s^3
}

kurtosis_hesapla <- function(x) {
  n <- length(x)
  ortalama <- mean(x)
  s <- sqrt(sum((x - ortalama)^2) / n)
  (sum((x - ortalama)^4) / n) / s^4 - 3
}

histogram_ciz <- function(x, baslik = NULL, x_ekseni_adi = "Puan") {
  ortalama <- mean(x)
  medyan <- median(x)
  mod_degeri <- mod_hesapla(x)

  cizgiler <- data.frame(
    tur = factor(c("Ortalama", "Medyan", "Mod"), levels = c("Ortalama", "Medyan", "Mod")),
    deger = c(ortalama, medyan, mod_degeri)
  )

  ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                    fill = "#5dade2", color = "white", alpha = 0.85) +
    geom_density(color = "#2c3e50", linewidth = 0.7) +
    geom_vline(data = cizgiler, aes(xintercept = deger, color = tur),
               linetype = "dashed", linewidth = 1) +
    scale_color_manual(values = c("Ortalama" = "#c0392b",
                                   "Medyan" = "#27ae60",
                                   "Mod" = "#8e44ad")) +
    labs(x = x_ekseni_adi, y = "Yoğunluk", title = baslik, color = "") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
}

ESIK_R_MUKEMMEL <- 0.40
ESIK_R_IYI <- 0.30
ESIK_R_SINIRDA <- 0.20
ESIK_P_COK_ZOR <- 0.20
ESIK_P_COK_KOLAY <- 0.80

BASARI_ESIK_YUKSEK <- 0.70
BASARI_ESIK_ORTA <- 0.50

madde_analizi_hesapla <- function(matris_0_1) {
  matris_0_1 <- as.matrix(matris_0_1)
  k <- ncol(matris_0_1)
  toplam <- rowSums(matris_0_1, na.rm = TRUE)

  p_degerleri <- colMeans(matris_0_1, na.rm = TRUE)
  r_degerleri <- sapply(seq_len(k), function(j) {
    diger_toplam <- toplam - matris_0_1[, j]
    suppressWarnings(cor(matris_0_1[, j], diger_toplam, use = "pairwise.complete.obs"))
  })

  karar <- vapply(r_degerleri, function(r) {
    if (is.na(r)) "Hesaplanamadı (sabit madde)"
    else if (r >= ESIK_R_MUKEMMEL) "TUT (çok iyi ayırt edici)"
    else if (r >= ESIK_R_IYI) "TUT (iyi, küçük revizyon önerilebilir)"
    else if (r >= ESIK_R_SINIRDA) "REVİZE ET (sınırda ayırt edicilik)"
    else "AT (yetersiz/negatif ayırt edicilik)"
  }, character(1))

  not_uyari <- vapply(p_degerleri, function(p) {
    if (is.na(p)) ""
    else if (p < ESIK_P_COK_ZOR) "Çok zor - r değerini baskılamış olabilir"
    else if (p > ESIK_P_COK_KOLAY) "Çok kolay - r değerini baskılamış olabilir"
    else ""
  }, character(1))

  madde_isimleri <- colnames(matris_0_1)
  if (is.null(madde_isimleri)) madde_isimleri <- paste0("Madde", seq_len(k))

  data.frame(
    Madde = madde_isimleri,
    p = round(p_degerleri, 3),
    r = round(r_degerleri, 3),
    Karar = karar,
    Not = not_uyari,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

ikili_matris_mi <- function(matris) {
  degerler <- unique(as.vector(as.matrix(matris)))
  degerler <- degerler[!is.na(degerler)]
  all(degerler %in% c(0, 1))
}

sablon_veri_uret <- function(n, a_vektor, b_vektor, seed) {
  set.seed(seed)
  k <- length(a_vektor)
  theta <- rnorm(n)
  matris <- matrix(NA_real_, n, k)
  for (j in seq_len(k)) {
    olasilik <- plogis(a_vektor[j] * theta - b_vektor[j])
    matris[, j] <- rbinom(n, 1, olasilik)
  }
  colnames(matris) <- paste0("Madde", seq_len(k))
  as.data.frame(matris)
}

SABLON_IYI <- sablon_veri_uret(
  n = 200,
  a_vektor = c(1.4, 1.6, 1.8, 1.5, 1.3, 1.7, 1.6, 1.4, 1.9, 1.5),
  b_vektor = c(-1.0, -0.6, -0.2, 0.2, 0.6, -0.4, 0.0, 0.4, -0.8, 0.8),
  seed = 101
)

SABLON_SORUNLU <- sablon_veri_uret(
  n = 200,
  a_vektor = c(1.5, 0.0, 1.6, 1.4, 1.5, 1.3, 1.6, 1.4, 1.5, 1.6),
  b_vektor = c(-3.5, 0.0, -0.4, 0.3, -0.6, 0.5, -0.2, 0.6, -0.5, 0.2),
  seed = 202
)

SABLON_KARISIK <- sablon_veri_uret(
  n = 200,
  a_vektor = c(1.6, 1.5, 0.5, 1.4, 0.3, 1.7, 1.5, 0.6, 1.4, 1.6, 0.4, 1.5),
  b_vektor = c(-0.5, 0.2, 0.3, -0.3, 0.8, 0.0, -0.6, -0.2, 0.4, 2.2, 0.5, -0.1),
  seed = 303
)

tema <- bs_theme(bootswatch = "zephyr")

bilgi_kutusu <- function(..., ikon = "lightbulb", renk = "info") {
  div(
    class = paste0("border-start border-4 border-", renk, " bg-", renk, "-subtle",
                    " rounded-2 p-2 mb-2 small"),
    div(class = "d-flex align-items-start gap-2",
        icon(ikon, class = paste0("text-", renk, " mt-1")),
        div(...))
  )
}

deger_kutusu <- function(title, value, showcase = NULL, theme = NULL) {
  value_box(
    title = title, value = value,
    showcase = showcase, theme = theme,
    showcase_layout = showcase_left_center(width = 0.22, max_height = "26px"),
    class = "kompakt-kutu",
    height = "95px"
  )
}

kaydirac_ucdeger_etiketi <- function(sol, orta, sag) {
  div(class = "d-flex justify-content-between text-muted",
      style = "font-size: 0.75rem; margin-top: -0.5rem;",
      span(sol), span(orta), span(sag))
}

ornek_csv_butonu <- function(ek_sinif = "") {
  tags$a(
    href = "../ornek_madde_verisi.csv",
    download = "ornek_madde_verisi.csv",
    class = paste("btn btn-outline-secondary btn-sm", ek_sinif),
    icon("download"), " Örnek CSV İndir"
  )
}

genel_css <- tags$style(HTML("
  .kompakt-kutu .value-box-value { font-size: 1.3rem !important; }
  .kompakt-kutu .value-box-title { font-size: 0.75rem !important; margin-bottom: 0.15rem; }
  .kompakt-kutu .value-box-showcase i { font-size: 1.3rem !important; }
  .card { box-shadow: 0 1px 4px rgba(0,0,0,.07); border: 1px solid rgba(0,0,0,.05); }
  .card-header { font-weight: 600; font-size: .9rem; }
  pre.shiny-text-output {
    background: #f8f9fb; border: 1px solid #e9ecef; border-radius: .5rem;
    font-size: .8rem; padding: .6rem .8rem; white-space: pre-wrap;
  }
  table.table { font-size: .85rem; }
  table.table th { text-transform: none; }
  .bslib-sidebar-layout > .sidebar { background-color: #fbfcfe; }
  .kaydirac-sayisiz .irs-min, .kaydirac-sayisiz .irs-max,
  .kaydirac-sayisiz .irs-single { display: none; }
"))

ui <- page_navbar(
  title = "Eğitimde Ölçme ve Değerlendirme: İnteraktif Uygulama Modülü",
  theme = tema,
  fillable = FALSE,
  header = genel_css,

  nav_panel(
    title = "Veri Laboratuvarı",
    icon = icon("chart-column"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always",
        width = 340,
        card_header(class = "bg-primary-subtle", "Kontrol Paneli: Dağılım Parametreleri"),
        p(class = "text-muted small",
          "Kaydırıcıları kullanarak çarpıklık ve basıklık değerlerini değiştirebilir, test puanlarının teorik dağılımını inceleyebilirsiniz."),
        div(class = "kaydirac-sayisiz",
            sliderInput("m1_carpiklik", "Çarpıklık:",
                        min = -0.4333, max = 0.4333, value = 0, step = 0.0217, ticks = FALSE),
            kaydirac_ucdeger_etiketi("Sola çarpık", "Simetrik", "Sağa çarpık")),
        div(class = "text-muted mb-2", style = "font-size: 0.8rem;",
            "→ Bu örneklemde gerçekleşen çarpıklık: ", strong(textOutput("m1_carpiklik_gercek", inline = TRUE))),
        div(class = "kaydirac-sayisiz",
            sliderInput("m1_basiklik", "Basıklık:",
                        min = -10, max = 10, value = 0, step = 0.5, ticks = FALSE),
            kaydirac_ucdeger_etiketi("Basık", "Normal", "Sivri")),
        div(class = "text-muted mb-2", style = "font-size: 0.8rem;",
            "→ Bu örneklemde gerçekleşen basıklık: ", strong(textOutput("m1_basiklik_gercek", inline = TRUE))),
        bilgi_kutusu(
          ikon = "circle-info", renk = "primary",
          strong("İpucu: "), "Kaydırıcıyı sağa çektikçe dağılım sivrileşir (uç değerler artar),",
          "sola çektikçe basıklaşır (uçlar azalır, dağılım yayvanlaşır);",
          "sıfırdayken normale yakın bir dağılım görürsünüz."
        ),
        hr(),
        uiOutput("m1_puan_secici_ui"),
        bilgi_kutusu(
          ikon = "bullseye", renk = "secondary",
          "Seçtiğiniz puan histogram üzerinde işaretlenir ve standart puan karşılıkları hesaplanır."
        )
      ),
      bilgi_kutusu(
        ikon = "graduation-cap", renk = "success",
        strong("Hoş geldiniz! "),
        "Bu uygulama, Ölçme ve Değerlendirme dersinin temel kavramlarını canlı olarak",
        "keşfetmeniz için hazırlandı. Sekmeleri soldan sağa sırayla izleyebilir ya da",
        "istediğiniz konuya doğrudan geçebilirsiniz. Her sekmede kaydırıcıları ve",
        "seçimleri değiştirerek sonuçların anında nasıl değiştiğini gözlemleyin."
      ),
      card(
        full_screen = TRUE,
        card_header("Canlı Histogram"),
        plotOutput("m1_histogram_plot", height = "400px")
      ),
      layout_columns(
        deger_kutusu("Mod", textOutput("m1_mod_deger"), icon("mountain"), "purple"),
        deger_kutusu("Medyan", textOutput("m1_medyan_deger"), icon("arrows-left-right"), "info"),
        deger_kutusu("Ortalama", textOutput("m1_ortalama_deger"), icon("scale-balanced"), "primary")
      ),
      bilgi_kutusu(
        ikon = "wand-magic-sparkles", renk = "warning",
        strong("Sıralama: "), textOutput("m1_siralama_metni", inline = TRUE),
        br(),
        textOutput("m1_yorum_metni"),
        br(),
        textOutput("m1_skew_kurt_metni")
      ),
      card(
        card_header("z / T / Yüzdelik Sıra Cetveli"),
        layout_columns(
          deger_kutusu("z puanı", textOutput("m1_z_deger"), theme = "secondary"),
          deger_kutusu("T puanı", textOutput("m1_t_deger"), theme = "secondary"),
          deger_kutusu("Yüzdelik Sıra", textOutput("m1_yuzdelik_deger"), theme = "secondary")
        )
      ),
      bilgi_kutusu(
        ikon = "flask", renk = "light",
        "Yüzdelik sıra, puanların yüzde kaçının seçilen puanın altında kaldığını",
        "gösterir (örn. 75 ise puanların %75'i bu puanın altındadır) ve bu sekmede",
        "normallik varsaymadan, doğrudan gözlenen veriden hesaplanır. Normallik",
        "varsayımına dayalı hesaplama için 'Standart Puanlar' sekmesini kullanabilirsiniz."
      )
    )
  ),

  nav_panel(
    title = "Madde Analizi",
    icon = icon("magnifying-glass-chart"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always",
        width = 340,
        card_header(class = "bg-primary-subtle", "Veri Kaynağı"),
        radioButtons("m2_kaynak", "Çalışma verisini seçin:",
                     choices = c("Hazır Şablon Kullan" = "sablon",
                                 "Kendi CSV Dosyamı Yükle" = "yukle"),
                     selected = "sablon"),
        conditionalPanel(
          condition = "input.m2_kaynak == 'sablon'",
          selectInput("m2_sablon_sec", "Şablon test türü:",
                      choices = c("İyi Ayırt Edici Test" = "iyi",
                                  "Sorunlu Test (Kusurlu Maddeler)" = "sorunlu",
                                  "Karma Nitelikli Test" = "karisik"))
        ),
        conditionalPanel(
          condition = "input.m2_kaynak == 'yukle'",
          bilgi_kutusu(
            ikon = "file-csv", renk = "secondary",
            "Format: Satırlar öğrenciyi, sütunlar maddeleri temsil etmelidir.",
            "Hücreler yalnızca 0 ve 1 değerlerini içermelidir (ID sütunu bulunmamalıdır)."
          ),
          fileInput("m2_veri_dosya", "CSV Dosyası Yükle:", accept = ".csv")
        ),
        ornek_csv_butonu("w-100 mt-2"),
        hr(),
        uiOutput("m2_madde_secici"),
        hr(),
        bilgi_kutusu(
          ikon = "scale-unbalanced", renk = "info",
          strong("Madde Ayırt Edicilik (r) Karar Eşikleri:"), br(),
          "r ≥ .40 → Tut (Mükemmel)", br(),
          ".30 ≤ r < .40 → Tut (İyi)", br(),
          ".20 ≤ r < .30 → Revize Et (Sınırda)", br(),
          "r < .20 → At (Yetersiz)", br(), br(),
          strong("Not: "), "p, maddeyi doğru cevaplayanların oranıdır;",
          "p büyüdükçe madde kolaylaşır."
        )
      ),
      layout_columns(
        deger_kutusu("Madde Güçlük İndeksi (p)", textOutput("m2_p_deger"), icon("gauge"), "primary"),
        deger_kutusu("Madde Ayırt Edicilik İndeksi (r)", textOutput("m2_r_deger"), icon("crosshairs"), "info")
      ),
      card(
        card_header("Karar: Madde Revizyon Durumu"),
        card_body(uiOutput("m2_karar_paneli"))
      ),
      card(
        card_header("Tüm Maddelerin Özet Tablosu"),
        p(class = "text-muted small",
          "p: Madde güçlüğü (doğru cevaplama oranı). r: Düzeltilmiş madde-toplam korelasyonu",
          "(parça-bütün yanlılığı giderilmiş ayırt edicilik indeksi)."),
        tableOutput("m2_ozet_tablo")
      )
    )
  ),

  nav_panel(
    title = "Güvenirlik ve Ölçme Hatası",
    icon = icon("shield-halved"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always",
        width = 340,
        card_header(class = "bg-primary-subtle", "Simülasyon Parametreleri"),
        sliderInput("n_madde", "Madde Sayısı (k):",
                    min = 2, max = 100, value = 20, step = 1),
        sliderInput("ort_r", "Maddeler Arası Ortalama Korelasyon (r̄):",
                    min = 0.05, max = 0.90, value = 0.30, step = 0.01),
        sliderInput("sd_test", "Test Toplam Puan Standart Sapması (SS):",
                    min = 1, max = 30, value = 15, step = 1),
        bilgi_kutusu(
          ikon = "square-root-variable", renk = "secondary",
          "Güvenirlik katsayısı (alpha), Spearman-Brown prensibine dayalı standardize",
          "formül ile hesaplanmaktadır:",
          withMathJax("$$\\alpha = \\frac{k \\bar{r}}{1+(k-1)\\bar{r}}$$"),
          withMathJax("$$\\text{ÖSH} = SS\\sqrt{1-\\alpha}$$")
        ),
        card(
          class = "border-warning border-3 shadow-sm",
          card_header(class = "bg-warning-subtle fw-bold",
                      icon("lightbulb"), " Aktif Öğrenme: Önce Tahmin Edin"),
          card_body(
            p(class = "small",
              "Mevcut ÖSH'nin yarıya inmesi için testte kaç madde olması",
              "gerektiğini tahmin edin (r̄ ve SS sabit tutularak)."),
            numericInput("tahmin_k", "Tahmin Edilen Madde Sayısı:",
                         value = 40, min = 2, max = 500, step = 1),
            actionButton("kontrol_et", "Tahmini Kontrol Et",
                         class = "btn-warning w-100 fw-bold"),
            br(), br(),
            verbatimTextOutput("tahmin_sonuc")
          )
        )
      ),
      layout_columns(
        deger_kutusu("Güncel Güvenirlik (α)", textOutput("alpha_deger"), icon("shield-heart"), "primary"),
        deger_kutusu("Ölçmenin Standart Hatası (ÖSH)", textOutput("sem_deger"), icon("ruler"), "danger")
      ),
      card(
        full_screen = TRUE,
        card_header("Madde Sayısı ve Hata İlişkisi"),
        plotOutput("guvenirlik_plot", height = "400px")
      ),
      bilgi_kutusu(
        ikon = "chart-line", renk = "warning",
        strong("Yorum: "),
        "Madde sayısı (k) arttıkça güvenirlik (alpha) katsayısı artar ve ölçmenin",
        "standart hatası (ÖSH) azalır. Ancak bu ilişki doğrusal değildir; belirli bir",
        "noktadan sonra fazladan eklenen maddelerin sağladığı psikometrik kazanç azalır."
      ),
      card(
        card_header("Gözlenen Veri Üzerinden Hesaplama (Ham Alpha)"),
        p(class = "text-muted small",
          "CSV formatındaki ham öğrenci yanıt matrisinizi yükleyerek doğrudan",
          "veriden hesaplanan güvenirlik katsayısını inceleyebilirsiniz.",
          "Elinizde veri yoksa örnek dosyayı indirip deneyebilirsiniz."),
        layout_columns(
          col_widths = c(8, 4),
          fileInput("veri_dosya", "CSV Dosyası Yükle:", accept = ".csv"),
          ornek_csv_butonu("mt-4")
        ),
        verbatimTextOutput("gercek_veri_sonuc")
      )
    )
  ),

  nav_panel(
    title = "Standart Puanlar ve Güven Aralığı",
    icon = icon("ruler-combined"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always",
        width = 340,
        card_header(class = "bg-primary-subtle", "Puan Bilgileri"),
        numericInput("ham_puan", "Ham Puan (X):", value = 70, step = 1),
        numericInput("ort_puan", "Grup Ortalaması (X̄):", value = 60, step = 1),
        bilgi_kutusu(
          ikon = "link", renk = "secondary",
          "Test standart sapması (SS) ve ölçmenin standart hatası (ÖSH) değerleri",
          "'Güvenirlik ve Ölçme Hatası' sekmesinden dinamik olarak aktarılmaktadır."
        ),
        checkboxInput("gercek_veri_kullan",
                      "Yüklenen CSV verisinin ÖSH değerini kullan",
                      value = FALSE),
        verbatimTextOutput("aktarilan_degerler"),
        hr(),
        radioButtons("ga_seviye", "Güven Aralığı (GA) Düzeyi:",
                     choices = c("%68 (±1 ÖSH)" = 1,
                                 "%95 (±1.96 ÖSH)" = 1.96,
                                 "%99 (±2.58 ÖSH)" = 2.58),
                     selected = 1.96)
      ),
      layout_columns(
        deger_kutusu("z puanı", textOutput("z_deger"), icon("z"), "primary"),
        deger_kutusu("T puanı", textOutput("t_deger"), icon("t"), "info"),
        deger_kutusu("Yüzdelik Sıra", textOutput("yuzdelik_deger"), icon("percent"), "secondary")
      ),
      card(
        full_screen = TRUE,
        card_header("Normal Dağılım ve Bireysel Puan Konumu"),
        plotOutput("normal_egri_plot", height = "400px")
      ),
      bilgi_kutusu(
        ikon = "cloud", renk = "info",
        strong("Yorum: "),
        "Yüzdelik sıra, normal dağılım varsayımı altında puanların yüzde kaçının",
        "bu puanın altında kaldığını gösterir. Grafikteki mavi gölgeli bant ise,",
        "adayın gerçek puanının (true score) %",
        textOutput("ga_yuzde_inline", inline = TRUE),
        " olasılıkla yer aldığı aralığı temsil etmektedir. Güvenirlik sekmesinde",
        "yapılan değişiklikler (örneğin k sayısının artırılması), ölçüm kesinliğini",
        "doğrudan etkileyerek bu güven aralığının daralmasına yol açacaktır."
      )
    )
  ),

  nav_panel(
    title = "Otomatik Analiz Modülü",
    icon = icon("laptop-code"),
    layout_sidebar(
      sidebar = sidebar(
        open = "always",
        width = 340,
        card_header(class = "bg-primary-subtle", "Veri Aktarımı"),
        p(class = "text-muted small",
          "Öğrenci yanıt matrisinizi yükleyerek test ve madde bazlı tüm özet",
          "istatistikleri tek bir arayüzde otomatik olarak elde edebilirsiniz."),
        bilgi_kutusu(
          ikon = "circle-info", renk = "secondary",
          "İkili (0-1) puanlanan testler için madde istatistikleri ve KR-20",
          "hesaplanır. Çoklu kategorik puanlamalarda (Likert tipi) ise yalnızca",
          "test düzeyi istatistikler ve ham Cronbach Alpha raporlanır."
        ),
        fileInput("m4_veri_dosya", "Veri Matrisi (CSV) Yükle:", accept = ".csv"),
        ornek_csv_butonu("w-100 mb-2"),
        uiOutput("m4_durum_metni")
      ),
      uiOutput("m4_panel")
    )
  )
)

server <- function(input, output, session) {

  m1_veri <- reactive({
    set.seed(42)
    dagilim_uret(n = 5000, carpiklik = input$m1_carpiklik,
                 basiklik_duzey = input$m1_basiklik,
                 hedef_ortalama = 60, hedef_sd = 15)
  })

  m1_istatistikler <- reactive({
    x <- m1_veri()
    list(
      ortalama = mean(x),
      medyan = median(x),
      mod = mod_hesapla(x),
      carpiklik_ornek = skewness_hesapla(x),
      basiklik_ornek = kurtosis_hesapla(x)
    )
  })

  m1_siralama_yorumu <- reactive({
    ist <- m1_istatistikler()
    degerler <- c(Mod = ist$mod, Medyan = ist$medyan, Ortalama = ist$ortalama)
    sirali <- sort(degerler)
    isimler <- names(sirali)
    siralama_metni <- paste(sprintf("%s (%.1f)", isimler, sirali), collapse = " < ")

    yorum <- if (ist$carpiklik_ornek > 0.15) {
      "Sağa çarpık (pozitif çarpıklık): uzun kuyruk sağda, ortalama mod ve medyanın sağına çekiliyor."
    } else if (ist$carpiklik_ornek < -0.15) {
      "Sola çarpık (negatif çarpıklık): uzun kuyruk solda, ortalama mod ve medyanın soluna çekiliyor."
    } else {
      "Yaklaşık simetrik: mod, medyan ve ortalama birbirine yakın."
    }
    list(siralama = siralama_metni, yorum = yorum)
  })

  output$m1_mod_deger <- renderText(sprintf("%.2f", m1_istatistikler()$mod))
  output$m1_medyan_deger <- renderText(sprintf("%.2f", m1_istatistikler()$medyan))
  output$m1_ortalama_deger <- renderText(sprintf("%.2f", m1_istatistikler()$ortalama))

  output$m1_siralama_metni <- renderText(m1_siralama_yorumu()$siralama)
  output$m1_yorum_metni <- renderText(m1_siralama_yorumu()$yorum)

  output$m1_carpiklik_gercek <- renderText(sprintf("%.2f", m1_istatistikler()$carpiklik_ornek))
  output$m1_basiklik_gercek <- renderText(sprintf("%.2f", m1_istatistikler()$basiklik_ornek))

  output$m1_skew_kurt_metni <- renderText({
    ist <- m1_istatistikler()
    sprintf("Örneklem çarpıklığı = %.2f  |  Örneklem basıklığı = %.2f (normal dağılımda 0)",
            ist$carpiklik_ornek, ist$basiklik_ornek)
  })

  output$m1_puan_secici_ui <- renderUI({
    x <- m1_veri()
    alt <- max(0, floor(min(x)))
    ust <- ceiling(max(x))
    mevcut <- isolate(input$m1_secilen_puan)
    baslangic <- if (!is.null(mevcut) && is.finite(mevcut) &&
                       mevcut >= alt && mevcut <= ust) mevcut else round(mean(x))
    sliderInput("m1_secilen_puan", "Cetvel için bir ham puan seçin:",
                min = alt, max = ust, value = baslangic, step = 1)
  })

  m1_z_t_yuzdelik <- reactive({
    x <- m1_veri()
    req(input$m1_secilen_puan)
    puan <- input$m1_secilen_puan
    ort <- mean(x)
    s <- sd(x)
    z <- (puan - ort) / s
    list(
      z = z,
      t = 50 + 10 * z,
      yuzdelik = mean(x <= puan) * 100
    )
  })

  output$m1_z_deger <- renderText(sprintf("%.2f", m1_z_t_yuzdelik()$z))
  output$m1_t_deger <- renderText(sprintf("%.1f", m1_z_t_yuzdelik()$t))
  output$m1_yuzdelik_deger <- renderText(sprintf("%.1f", m1_z_t_yuzdelik()$yuzdelik))

  output$m1_histogram_plot <- renderPlot({
    x <- m1_veri()
    grafik <- histogram_ciz(x, x_ekseni_adi = "Puan")
    req(input$m1_secilen_puan)
    grafik +
      geom_vline(xintercept = input$m1_secilen_puan, color = "#f39c12", linewidth = 1.1) +
      annotate("text", x = input$m1_secilen_puan, y = Inf, vjust = 1.5,
               label = "Seçilen Puan", color = "#f39c12")
  })

  m2_yuklenen <- reactiveValues(veri = NULL, hata = NULL)

  observeEvent(input$m2_veri_dosya, {
    tryCatch({
      df <- read.csv(input$m2_veri_dosya$datapath)
      df_sayisal <- df[sapply(df, is.numeric)]

      if (ncol(df_sayisal) < 2) {
        m2_yuklenen$hata <- "En az 2 sayısal madde sütunu gereklidir."
        m2_yuklenen$veri <- NULL
        return()
      }
      if (nrow(df_sayisal) < 2) {
        m2_yuklenen$hata <- "En az 2 öğrenci (satır) gereklidir."
        m2_yuklenen$veri <- NULL
        return()
      }
      if (!ikili_matris_mi(df_sayisal)) {
        m2_yuklenen$hata <- "Madde analizi için tüm hücreler yalnızca 0 ya da 1 olmalıdır."
        m2_yuklenen$veri <- NULL
        return()
      }

      m2_yuklenen$veri <- df_sayisal
      m2_yuklenen$hata <- NULL

    }, error = function(e) {
      m2_yuklenen$hata <- paste("Dosya okunamadı veya hesaplama hatası:", e$message)
      m2_yuklenen$veri <- NULL
    })
  })

  m2_aktif_veri <- reactive({
    if (input$m2_kaynak == "sablon") {
      switch(input$m2_sablon_sec,
             iyi = SABLON_IYI,
             sorunlu = SABLON_SORUNLU,
             karisik = SABLON_KARISIK)
    } else {
      m2_yuklenen$veri
    }
  })

  output$m2_madde_secici <- renderUI({
    if (input$m2_kaynak == "yukle") {
      if (!is.null(m2_yuklenen$hata)) {
        return(div(class = "text-danger fw-semibold small", m2_yuklenen$hata))
      }
      if (is.null(m2_yuklenen$veri)) {
        return(helpText("Lütfen bir CSV dosyası yükleyin."))
      }
    }
    veri <- m2_aktif_veri()
    req(veri)
    selectInput("m2_madde_sec", "Analiz Edilecek Madde:", choices = colnames(veri))
  })

  m2_analiz_sonuc <- reactive({
    veri <- m2_aktif_veri()
    req(veri)
    madde_analizi_hesapla(veri)
  })

  m2_secili_satir <- reactive({
    req(input$m2_madde_sec)
    sonuc <- m2_analiz_sonuc()
    sonuc[sonuc$Madde == input$m2_madde_sec, ]
  })

  output$m2_p_deger <- renderText({
    satir <- m2_secili_satir()
    req(nrow(satir) == 1)
    sprintf("%.3f", satir$p)
  })

  output$m2_r_deger <- renderText({
    satir <- m2_secili_satir()
    req(nrow(satir) == 1)
    sprintf("%.3f", satir$r)
  })

  output$m2_karar_paneli <- renderUI({
    satir <- m2_secili_satir()
    req(nrow(satir) == 1)

    alert_sinifi <- if (grepl("^TUT", satir$Karar)) "alert-success"
                    else if (grepl("^REVİZE", satir$Karar)) "alert-warning"
                    else "alert-danger"
    alert_ikonu <- if (grepl("^TUT", satir$Karar)) "circle-check"
                   else if (grepl("^REVİZE", satir$Karar)) "triangle-exclamation"
                   else "circle-xmark"

    div(
      class = paste("alert mb-0", alert_sinifi),
      h5(class = "alert-heading mb-0", icon(alert_ikonu), " ", satir$Karar),
      if (nchar(satir$Not) > 0)
        p(class = "small mt-2 mb-0", strong("Not: "), satir$Not)
    )
  })

  output$m2_ozet_tablo <- renderTable({
    sonuc <- m2_analiz_sonuc()
    sonuc[, setdiff(names(sonuc), "Not")]
  }, striped = TRUE, bordered = TRUE)

  alpha_reactive <- reactive({
    k <- input$n_madde
    r <- input$ort_r
    (k * r) / (1 + (k - 1) * r)
  })

  sem_reactive <- reactive({
    input$sd_test * sqrt(1 - alpha_reactive())
  })

  output$alpha_deger <- renderText({
    sprintf("%.3f", alpha_reactive())
  })

  output$sem_deger <- renderText({
    sprintf("%.2f", sem_reactive())
  })

  output$guvenirlik_plot <- renderPlot({
    k_seq <- 2:100
    r <- input$ort_r
    sd_test <- input$sd_test

    alpha_seq <- (k_seq * r) / (1 + (k_seq - 1) * r)
    sem_seq <- sd_test * sqrt(1 - alpha_seq)

    df <- data.frame(k = k_seq, alpha = alpha_seq, sem = sem_seq)

    olcek <- max(df$sem) / max(df$alpha)

    ggplot(df, aes(x = k)) +
      geom_line(aes(y = alpha, color = "Güvenirlik (α)"), linewidth = 1.1) +
      geom_line(aes(y = sem / olcek, color = "Ölçme Hatası (ÖSH)"), linewidth = 1.1) +
      geom_vline(xintercept = input$n_madde, linetype = "dashed", color = "gray40") +
      scale_y_continuous(
        name = "Güvenirlik (α)",
        sec.axis = sec_axis(~ . * olcek, name = "ÖSH")
      ) +
      scale_color_manual(values = c("Güvenirlik (α)" = "#2c3e50",
                                     "Ölçme Hatası (ÖSH)" = "#c0392b")) +
      labs(x = "Madde Sayısı (k)", color = "") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
  })

  tahmin_metni <- eventReactive(input$kontrol_et, {
    validate(need(!is.null(input$tahmin_k) && is.finite(input$tahmin_k),
                  "Lütfen önce bir madde sayısı tahmini girin."))
    r <- input$ort_r
    alpha_mevcut <- alpha_reactive()
    tahmin_k_o_an <- round(input$tahmin_k)

    alpha_hedef <- 1 - 0.25 * (1 - alpha_mevcut)

    k_gerekli <- (alpha_hedef * (1 - r)) / (r * (1 - alpha_hedef))
    k_gerekli <- ceiling(k_gerekli)

    fark <- tahmin_k_o_an - k_gerekli

    if (alpha_hedef >= 0.999) {
      paste0("Bu r̄ değeriyle ÖSH'yi yarıya indirmek pratik olarak",
             " çok yüksek madde sayısı gerektiriyor (k > 500).",
             " Daha yüksek bir r̄ ile tekrar deneyin.")
    } else {
      sprintf(paste0("Gerekli madde sayısı: k ≈ %d\n",
                     "Sizin tahmininiz: %d\n",
                     "Fark: %+d madde\n\n",
                     "%s"),
              k_gerekli, tahmin_k_o_an, fark,
              if (abs(fark) <= 5) "Tahmininiz oldukça isabetli." else
                "Formülün doğrusal olmadığını unutmayın; azalan verim ilkesi geçerli.")
    }
  })

  output$tahmin_sonuc <- renderText(tahmin_metni())

  gercek_veri <- reactiveValues(alpha = NULL, sem = NULL, sd = NULL, k = NULL, hata = NULL)

  observeEvent(input$veri_dosya, {
    tryCatch({
      df <- read.csv(input$veri_dosya$datapath)

      df_sayisal <- df[sapply(df, is.numeric)]

      if (ncol(df_sayisal) < 2) {
        gercek_veri$hata <- "En az 2 sayısal madde sütunu gereklidir."
        gercek_veri$alpha <- NULL
        return()
      }
      if (nrow(df_sayisal) < 2) {
        gercek_veri$hata <- "En az 2 öğrenci (satır) gereklidir."
        gercek_veri$alpha <- NULL
        return()
      }

      alpha_deger <- ham_alpha_hesapla(df_sayisal)

      toplam_puan <- rowSums(df_sayisal, na.rm = TRUE)
      sd_gercek <- sd(toplam_puan, na.rm = TRUE)
      sem_gercek <- sd_gercek * sqrt(1 - alpha_deger)

      gercek_veri$alpha <- alpha_deger
      gercek_veri$sd <- sd_gercek
      gercek_veri$sem <- sem_gercek
      gercek_veri$k <- ncol(df_sayisal)
      gercek_veri$hata <- NULL

    }, error = function(e) {
      gercek_veri$hata <- paste("Dosya okunamadı veya hesaplama hatası:", e$message)
      gercek_veri$alpha <- NULL
    })
  })

  output$gercek_veri_sonuc <- renderText({
    if (is.null(input$veri_dosya)) {
      return("Henüz dosya yüklenmedi. Simülasyon değerleri kullanılıyor.")
    }
    if (!is.null(gercek_veri$hata)) {
      return(gercek_veri$hata)
    }
    sprintf(paste0("Madde sayısı (gerçek veri) = %d\n",
                   "Cronbach alpha - ham (gerçek veri) = %.3f\n",
                   "Toplam puan SS (gerçek veri) = %.2f\n",
                   "ÖSH (gerçek veri) = %.2f\n\n",
                   "Not: Bu değer, kovaryans matrisinden hesaplanan HAM",
                   "(raw) alpha'dır; yukarıdaki simülasyon eğrisi ise",
                   "STANDARDİZE alpha formülüne dayanır. İki tanım genellikle",
                   "birbirine yakın ama özdeş değildir."),
            gercek_veri$k, gercek_veri$alpha, gercek_veri$sd, gercek_veri$sem)
  })

  aktif_sem <- reactive({
    if (isTRUE(input$gercek_veri_kullan) && !is.null(gercek_veri$sem)) {
      gercek_veri$sem
    } else {
      sem_reactive()
    }
  })

  aktif_sd <- reactive({
    if (isTRUE(input$gercek_veri_kullan) && !is.null(gercek_veri$sd)) {
      gercek_veri$sd
    } else {
      input$sd_test
    }
  })

  output$aktarilan_degerler <- renderText({
    kaynak <- if (isTRUE(input$gercek_veri_kullan) && !is.null(gercek_veri$sem)) {
      "Gerçek Veri (Güvenirlik sekmesindeki CSV)"
    } else {
      "Simülasyon (Güvenirlik sekmesi)"
    }
    sprintf("Kaynak: %s\nSS  = %.2f\nÖSH = %.2f",
            kaynak, aktif_sd(), aktif_sem())
  })

  z_reactive <- reactive({
    req(input$ham_puan, input$ort_puan)
    (input$ham_puan - input$ort_puan) / aktif_sd()
  })

  output$z_deger <- renderText({
    sprintf("%.2f", z_reactive())
  })

  output$t_deger <- renderText({
    sprintf("%.1f", 50 + 10 * z_reactive())
  })

  output$yuzdelik_deger <- renderText({
    sprintf("%.1f", pmin(99.9, pmax(0.1, pnorm(z_reactive()) * 100)))
  })

  output$ga_yuzde_inline <- renderText({
    switch(input$ga_seviye,
           "1" = "68",
           "1.96" = "95",
           "2.58" = "99")
  })

  output$normal_egri_plot <- renderPlot({
    req(input$ham_puan, input$ort_puan)
    sd_kullanilan <- aktif_sd()
    sem <- aktif_sem()
    katsayi <- as.numeric(input$ga_seviye)
    alt_sinir <- input$ham_puan - katsayi * sem
    ust_sinir  <- input$ham_puan + katsayi * sem

    min_x <- min(input$ort_puan - 4 * sd_kullanilan, alt_sinir - 5)
    max_x <- max(input$ort_puan + 4 * sd_kullanilan, ust_sinir + 5)

    x_ekseni <- seq(min_x, max_x, length.out = 500)
    yogunluk <- dnorm(x_ekseni, mean = input$ort_puan, sd = sd_kullanilan)
    df <- data.frame(x = x_ekseni, y = yogunluk)

    ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1) +
      geom_ribbon(data = subset(df, x >= alt_sinir & x <= ust_sinir),
                  aes(ymin = 0, ymax = y), fill = "#3498db", alpha = 0.35) +
      geom_vline(xintercept = input$ham_puan, color = "#c0392b",
                 linetype = "dashed", linewidth = 1) +
      annotate("text", x = input$ham_puan, y = max(df$y) * 1.05,
               label = "Gözlenen Puan", color = "#c0392b") +
      labs(x = "Ham Puan", y = "Yoğunluk") +
      theme_minimal(base_size = 14)
  })

  m4_yuklenen <- reactiveValues(veri = NULL, hata = NULL)

  observeEvent(input$m4_veri_dosya, {
    tryCatch({
      df <- read.csv(input$m4_veri_dosya$datapath)
      df_sayisal <- df[sapply(df, is.numeric)]

      if (ncol(df_sayisal) < 2) {
        m4_yuklenen$hata <- "En az 2 sayısal madde sütunu gereklidir."
        m4_yuklenen$veri <- NULL
        return()
      }
      if (nrow(df_sayisal) < 2) {
        m4_yuklenen$hata <- "En az 2 öğrenci (satır) gereklidir."
        m4_yuklenen$veri <- NULL
        return()
      }

      m4_yuklenen$veri <- df_sayisal
      m4_yuklenen$hata <- NULL

    }, error = function(e) {
      m4_yuklenen$hata <- paste("Dosya okunamadı veya hesaplama hatası:", e$message)
      m4_yuklenen$veri <- NULL
    })
  })

  m4_ozet <- reactive({
    veri <- m4_yuklenen$veri
    req(veri)
    toplam <- rowSums(veri, na.rm = TRUE)
    ikili <- ikili_matris_mi(veri)
    list(
      veri = veri,
      toplam = toplam,
      ikili = ikili,
      ortalama = mean(toplam),
      sd = sd(toplam),
      minimum = min(toplam),
      maksimum = max(toplam),
      alpha = ham_alpha_hesapla(veri),
      etiket = if (ikili) "KR-20" else "Cronbach Alpha (ham)",
      k = ncol(veri)
    )
  })

  output$m4_durum_metni <- renderUI({
    if (is.null(input$m4_veri_dosya)) {
      return(helpText("Henüz dosya yüklenmedi."))
    }
    if (!is.null(m4_yuklenen$hata)) {
      return(div(class = "text-danger fw-semibold small", m4_yuklenen$hata))
    }
    ozet <- m4_ozet()
    div(class = "text-success fw-semibold small",
        icon("circle-check"), " ",
        sprintf("%d öğrenci, %d madde başarıyla yüklendi.", nrow(ozet$veri), ozet$k))
  })

  output$m4_panel <- renderUI({
    req(m4_yuklenen$veri)
    ozet <- m4_ozet()

    basari_kutusu <- NULL
    if (ozet$ikili) {
      basari_yuzdesi <- ozet$ortalama / ozet$k * 100
      basari_yorumu <- if (basari_yuzdesi >= BASARI_ESIK_YUKSEK * 100) {
        "YÜKSEK"
      } else if (basari_yuzdesi >= BASARI_ESIK_ORTA * 100) {
        "ORTA"
      } else {
        "DÜŞÜK"
      }
      basari_renk <- switch(basari_yorumu,
                            "YÜKSEK" = "success",
                            "ORTA" = "warning",
                            "DÜŞÜK" = "danger")
      basari_kutusu <- bilgi_kutusu(
        ikon = "graduation-cap", renk = basari_renk,
        strong("Sınıfın Başarı Düzeyi: "), basari_yorumu,
        sprintf(" (ortalama başarı yüzdesi = %%%.1f)", basari_yuzdesi),
        br(),
        em("Bu basit, öğretim amaçlı bir sınıflandırmadır (>%70 yüksek,",
           "%50-70 orta, <%50 düşük); psikometrik bir standart değildir.",
           "Sınıf başarısı sadece ortalamaya değil; testin güçlük düzeyine ve",
           "puanların dağılımına bağlı olarak da yorumlanmalıdır.")
      )
    }

    tagList(
      layout_columns(
        deger_kutusu("Test Ortalaması", sprintf("%.2f", ozet$ortalama),
                     icon("scale-balanced"), "primary"),
        deger_kutusu("Standart Sapma (SS)", sprintf("%.2f", ozet$sd),
                     icon("arrows-left-right"), "info"),
        deger_kutusu("Min - Maks", sprintf("%.0f - %.0f", ozet$minimum, ozet$maksimum),
                     icon("ruler-horizontal"), "secondary"),
        deger_kutusu(ozet$etiket, sprintf("%.3f", ozet$alpha),
                     icon("shield-heart"), "purple")
      ),
      basari_kutusu,
      card(
        full_screen = TRUE,
        card_header("Puan Dağılımı"),
        plotOutput("m4_histogram_plot", height = "350px")
      ),
      if (ozet$ikili) {
        card(
          card_header("Madde Bazlı Analiz (p / r / Karar)"),
          tableOutput("m4_madde_tablosu")
        )
      } else {
        bilgi_kutusu(
          ikon = "circle-info", renk = "secondary",
          "Not: Madde bazlı p/r analiz tablosu yalnızca 0/1 (ikili)",
          "puanlanan maddeler için hesaplanır. Yüklenen veri ikili",
          "olmadığı için bu tablo gösterilmiyor."
        )
      }
    )
  })

  output$m4_histogram_plot <- renderPlot({
    ozet <- m4_ozet()
    histogram_ciz(ozet$toplam, x_ekseni_adi = "Toplam Puan")
  })

  output$m4_madde_tablosu <- renderTable({
    ozet <- m4_ozet()
    req(ozet$ikili)
    madde_analizi_hesapla(ozet$veri)
  }, striped = TRUE, bordered = TRUE)
}

shinyApp(ui = ui, server = server)
