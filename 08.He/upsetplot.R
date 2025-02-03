library(ComplexUpset)
data<-read.csv("Orthogroups.frame.csv",header = T)
upset(data,intersect = colnames(data)[2:8],
      set_sizes = F,
      sort_intersections_by=c('degree', 'cardinality'),
      base_annotations=list (
          'Intersection size'=intersection_size(
              text = list(
                  vjust=0.5,
                  hjust=-0.1,
                  angle=90)
              )
          ),
      queries = list(upset_query(intersect = c("Maximus","Laperouse","Vlamingh","Buff","Yeti","Stirling","Clipper"),
                                 color="red",fill="red",
                                 only_components = c('intersections_matrix', 'Intersection size')),
                     upset_query(intersect = "Maximus",color="#91CF60",fill="#91CF60",only_components = c('intersections_matrix', 'Intersection size')),
                     upset_query(intersect = "Laperouse",color="#91CF60",fill="#91CF60",only_components = c('intersections_matrix', 'Intersection size')),
                     upset_query(intersect = "Vlamingh",color="#91CF60",fill="#91CF60",only_components = c('intersections_matrix', 'Intersection size')),
                     upset_query(intersect = "Buff",color="#91CF60",fill="#91CF60",only_components = c('intersections_matrix', 'Intersection size')),
                     upset_query(intersect = "Yeti",color="#91CF60",fill="#91CF60",only_components = c('intersections_matrix', 'Intersection size')),
                     upset_query(intersect = "Stirling",color="#91CF60",fill="#91CF60",only_components = c('intersections_matrix', 'Intersection size')),
                     upset_query(intersect = "Clipper",color="#91CF60",fill="#91CF60",only_components = c('intersections_matrix', 'Intersection size'))
                     )
      )
