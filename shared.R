# ks - kernel smoothing
# všechny varianty ks navíc mají ještě "null" a podvarianty dupl+uniq
versionNames <- function() {
    vn <- list(
        "1" = "TGOB_5000_kss_CZ", "2" = "TGOB_5000_kss_SSOSOBclip",
        "3" = "SSOSOB_5000_kss_CZ", "4" = "SSOSOB_SSOSOB_kss_CZ", "5" = "SSOSOB_P_kss_CZ",
        "6" = "SSOSOB_5000_kss_SSOSOBclip", "7" = "SSOSOB_SSOSOB_kss_SSOSOBclip", "8" = "SSOSOB_P_kss_SSOSOBclip",
        "9" = "TGOB_TGOB_fix_x", "10" = "SSOSOB_SSOSOB_fix_x",
        "11" = "TGOB_5000_rss_x", "12" = "TGOB_SSOSOB_rss_x", "13" = "TGOB_P_rss_x",
        "14" = "SSOSOB_5000_rss_x", "15" = "SSOSOB_SSOSOB_rss_x", "17" = "SSOSOB_P_rss_x",
        "17" = "SSOSOB_5000_rss_x_thin", "18" = "SSOSOB_SSOSOB_rss_x_thin", "19" = "SSOSOB_P_rss_x_thin",
        "20" = "TGOB_5000_rss_x_thin", "21" = "TGOB_SSOSOB_rss_x_thin", "22" = "TGOB_P_rss_x_thin"
    )
    return(vn)
}