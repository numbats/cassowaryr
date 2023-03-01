library(hexSticker)
library(showtext)
font_add_google("Love Ya Like A Sister", "gochi")
showtext_auto()

imgurl <- "heximage.png"

s <- sticker(imgurl, package="CassowaryR", p_color = "firebrick", p_size=28, p_x=1,
             p_y = 1.6, h_fill="white",h_color="darkblue", p_family = "gochi",
             s_x=1.1, s_y=0.9, s_width = 0.9, filename="hexfile.png", dpi = 600)

plot(s)

usethis::use_logo("hexfile.png")

# 2nd choice
#font_add_google("Monofett", "gochi")
#showtext_auto()
#imgurl <- system.file("figures/heximage2c.png", package="hexSticker")
#s <- sticker(imgurl, package="CassowaryR", p_color = "darkblue", p_size=4.5, p_x=1.4,
#             p_y = 0.3, h_fill="white",h_color="darkblue", p_family = "gochi",
#             s_x=1.1, s_y=1, s_width = 0.95, filename="hexfile.png", angle=30)
#plot(s)
