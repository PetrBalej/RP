# ks - kernel smoothing
# všechny varianty ks navíc mají ještě "null" a podvarianty dupl+uniq
versionNames <- function() {
    vn <- list(
        "1" = "TGOB_ks_CZ_5000", "2" = "TGOB_ks_SSOOGclip_5000",
        "3" = "SSOOG_ks_CZ_5000", "4" = "SSOOG_ks_CZ_SSOOG", "5" = "SSOOG_ks_CZ_P",
        "6" = "SSOOG_ks_SSOOGclip_5000", "7" = "SSOOG_ks_SSOOGclip_SSOOG", "8" = "SSOOG_ks_SSOOGclip_P",
        "9" = "TGOB_x_CZ_P+TGOB_x_SSOOGclip_P" # dodatečně rozdělit?
    )
    return(vn)
}