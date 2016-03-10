downtriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/100 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size*0.35, vertex.size, vertex.size*0.7,vertex.size,vertex.size*0.35,vertex.size*0.34,vertex.size*0.7,vertex.size*0.34),
          add=TRUE, inches=FALSE)
}

uptriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/100 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size*0.35, vertex.size*0.34, vertex.size*0.7,vertex.size*0.34,vertex.size*0.35,vertex.size,vertex.size*0.7,vertex.size),
          add=TRUE, inches=FALSE)
}


