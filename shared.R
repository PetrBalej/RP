# ks - kernel smoothing
# všechny varianty ks navíc mají ještě "null" a podvarianty dupl+uniq
versionNames <- function() {
    vn <- list(
        "1" = list("short" = "TGOB_5000_kss_CZ_x", "desc" = "tradiční TGOB, jako BG generování smoothing kernelem 5000 bodů", "group" = "1"),
        "2" = list("short" = "TGOB_5000_kss_SSOSOBclip_x", "desc" = "tradiční TGOB, jako BG generování smoothing kernelem 5000 bodů do oblasti SSOSOB", "group" = "1i"),
        ###
        "3" = list("short" = "SSOSOB_5000_kss_CZ_x", "desc" = "SSOSBG modifikace 1, jako BG generování smoothing kernelem 5000 bodů", "group" = "1I"),
        "4" = list("short" = "SSOSOB_SSOSOB_kss_CZ_x", "desc" = "SSOSBG modifikace 1, jako BG generování smoothing kernelem (n=SSOSBG) bodů", "group" = "1I"),
        "5" = list("short" = "SSOSOB_P_kss_CZ_x", "desc" = "SSOSBG modifikace 1, jako BG generování smoothing kernelem (n=P, všechny unik. presence druhu) bodů", "group" = "1I"),
        ###
        "6" = list("short" = "SSOSOB_5000_kss_SSOSOBclip_x", "desc" = "SSOSBG modifikace 1+3, jako BG generování smoothing kernelem 5000 bodů do oblasti SSOSOB", "group" = "1I"),
        "7" = list("short" = "SSOSOB_SSOSOB_kss_SSOSOBclip_x", "desc" = "SSOSBG modifikace 1+4, jako BG generování smoothing kernelem (n=SSOSBG) bodů do oblasti SSOSOB", "group" = "1I"),
        "8" = list("short" = "SSOSOB_P_kss_SSOSOBclip_x", "desc" = "SSOSBG modifikace 1+5, jako BG generování smoothing kernelem (n=P, všechny unik. presence druhu) bodů do oblasti SSOSOB", "group" = "1I"),
        ###
        "9" = list("short" = "TGOB_TGOB_fix_x_x", "desc" = "tradiční TGOB, jako BG všechny unikátní pixely (~sites) presencí", "group" = "0"),
        "10" = list("short" = "SSOSOB_SSOSOB_fix_x_x", "desc" = "SSOSBG modifikace 9, BG jsou unikátní pixely SSOSBG", "group" = "0I"), # vzdálená podobnost generování BG kolem presencí, zde ale generiju BG do prozkoumaných oblastí, často úplně mimo presence druhu
        ###
        "11" = list("short" = "TGOB_5000_rss_x_x", "desc" = "tradiční TGOB, jako BG random subsampling 5000 presencí z TGOB", "group" = "0"),
        "12" = list("short" = "TGOB_SSOSOB_rss_x_x", "desc" = "SSOSBG modifikace 11, jako BG random subsampling presencí (n=SSOSBG) z TGOB", "group" = "0i"),
        "13" = list("short" = "TGOB_P_rss_x_x", "desc" = "P modifikace 11, jako BG random subsampling presencí (n=P, všechny unik. presence druhu) z TGOB", "group" = "0"),
        ###
        "14" = list("short" = "SSOSOB_5000_rss_x_x", "desc" = "SSOSBG modifikace 9, jako BG random subsampling 5000 presencí z SSOSBG", "group" = "0I"),
        "15" = list("short" = "SSOSOB_SSOSOB_rss_x_x", "desc" = "SSOSBG modifikace 9, jako BG random subsampling presencí (n=SSOSBG) z SSOSBG", "group" = "0I"),
        "16" = list("short" = "SSOSOB_P_rss_x_x", "desc" = "SSOSBG modifikace 9, jako BG random subsampling presencí (n=P, všechny unik. presence druhu) z SSOSBG", "group" = "0I"),
        ###
        "17" = list("short" = "SSOSOB_5000_rss_x_thin", "desc" = "totéž co 14, navíc thinning presencí", "group" = "0I"),
        "18" = list("short" = "SSOSOB_SSOSOB_rss_x_thin", "desc" = "totéž co 15, navíc thinning presencí", "group" = "0I"),
        "19" = list("short" = "SSOSOB_P_rss_x_thin", "desc" = "totéž co 16, navíc thinning presencí", "group" = "0I"),
        ###
        "20" = list("short" = "TGOB_5000_rss_x_thin", "desc" = "totéž co 11, navíc thinning presencí", "group" = "0"),
        "21" = list("short" = "TGOB_SSOSOB_rss_x_thin", "desc" = "totéž co 12, navíc thinning presencí", "group" = "0i"),
        "22" = list("short" = "TGOB_P_rss_x_thin", "desc" = "totéž co 13, navíc thinning presencí", "group" = "0")
    )
    return(vn)
}

versionNamesDf <- function() {
    vn <- versionNames()
    df.colnames <- c("id", "bg_source", "bg_n", "bg_sampling_method", "bg_region", "presence_thinning", "description", "group")
    first <- TRUE
    for (id in names(vn)) {
        df.temp <- as_tibble(t(c(id, strsplit(vn[[id]][["short"]], "_", fixed = TRUE)[[1]], vn[[id]][["desc"]], vn[[id]][["group"]])))
        colnames(df.temp) <- df.colnames
        if (first) {
            first <- FALSE
            df <- df.temp
        } else {
            df %<>% add_row(df.temp)
        }
    }
    return(df)
}