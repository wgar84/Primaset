require(devtools)
require(roxygen2)
require(shapes)

### test

avahi.test <- prosimian.raw [[1]] [[2]] [, , 1, 3]

shapes3d(avahi.test[, ])

avahi.glued <- glueSkull(avahi.test [1:28, ], avahi.test [29:56, ])

rownames(avahi.test [1:28, ])

shapes3d(avahi.glued)

### symmetries

rownames(prosimian.raw [[268]] [[2]] [1:28, , 1, 1])

