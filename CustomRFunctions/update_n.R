update_n <- function (x, showCategory) 
{
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  return(n)
}





## 2-25-19
## me playing with this function

x <- ego
update_n(x, showCategory = 5)

x <- ego[1:4]
update_n(x, showCategory = 5)

x <- ego
update_n(x, showCategory = "poop")


# what it does: sets the number of GO terms to display to 5 or less