# a <- rep(seq(10), 10)/10
# b <- rep(seq(10), each=10)/10
# sc <- scree(a, b)
# mst <- gen_mst(sc$del, sc$weights)
# morethantwo <- which(igraph::degree(mst)>=2)
# adjs <- mst[morethantwo]
# points <- sc$del$x[]
# origin <- sc$del$x[angs[1],]
# vects <- t(t(points)-origin)
#
# as_adjacency_matrix(mst)
# as_adjacency_matrix(mst)
#
# plot(a,b)
# sc_grid(a,b)
# draw_mst(a, b) + theme_minimal() +
#   geom_point(data = tibble(x=c(0.01269430),
#                            y = c(0.05228626)),
#              aes(x,y), colour="red")
