## The function computes month, day, and year from Julian days. 
"caldat" = function (julian) 
{
  igreg = 2299161
  julian <- trunc(julian)
  jalpha <- julian * 0
  ja <- julian * 0
  im <- (julian >= igreg)
  if (sum(im) > 0) {
    jalpha[im] <- trunc(((julian - 1867216) - 0.25)/36524.25)
    ja[im] <- julian + 1 + jalpha - trunc(0.25 * jalpha)
  }
  im <- (julian < igreg)
  if (sum(im) > 0) 
    ja[im] <- julian[im]
  jb <- ja + 1524
  jc <- trunc(6680 + ((jb - 2439870) - 122.1)/365.25)
  jd <- 365 * jc + trunc(0.25 * jc)
  je <- trunc((jb - jd)/30.6001)
  id <- jb - jd - trunc(30.6001 * je)
  mm <- je - 1
  im <- (mm > 12)
  if (sum(im) > 0) 
    mm[im] <- mm[im] - 12
  iyyy <- jc - 4715
  im <- (mm > 2)
  if (sum(im) > 0) 
    iyyy[im] <- iyyy[im] - 1
  im <- (iyyy <= 0)
  if (sum(im) > 0) 
    iyyy <- iyyy - 1
  caldat <- list(month = mm, day = id, year = iyyy)
  invisible(caldat)
}

## The function computes Julian days from month, day, and year. 
"julday" <- function (mm, id, iyyy)
{
  igreg <- 588829
  mm <- trunc(mm)
  id <- trunc(id)
  iyyy <- trunc(iyyy)
  im <- (iyyy == 0)
  if (sum(im, na.rm = TRUE) > 0)
    return("There is no year zero!")
  if ((length(mm) != length(id)) | (length(mm) != length(iyyy)) |
      (length(iyyy) != length(id)))
    return("The vectors must have same length!")
  im <- (iyyy < 0)
  if (sum(im) > 0)
    iyyy[im] <- iyyy[im] + 1
  jy <- mm * 0
  jm <- mm * 0
  ja <- mm * 0
  im <- (mm > 2)
  if (sum(im) > 0) {
    jy[im] <- iyyy[im]
    jm[im] <- mm[im] + 1
  }
  im <- (mm <= 2)
  if (sum(im) > 0) {
    jy[im] <- iyyy[im] - 1
    jm[im] <- mm[im] + 13
  }
  jul <- trunc(365.25 * jy) + trunc(30.6001 * jm) + id + 1720995
  im <- (id + 31 * (mm + 12 * iyyy) >= igreg)
  if (sum(im) > 0) {
    ja[im] <- trunc(0.01 * jy)
    jul[im] <- jul + 2 - ja[im] + trunc(0.25 * ja[im])
  }
  julday <- jul
  invisible(julday)
} 