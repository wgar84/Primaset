## plotar W

hyp.nominal <- c()
for (i in 1:39)
  {
    hyp.nominal [i] <-
      paste (colnames (Aux $ def.hyp) [1:6] [which (Aux $ def.hyp [i, 1:6] == 1)],
             collapse = '/')
  }

hyp.nominal [1] <- 'None'

hyp.nominal <- hyp.nominal[-1]

hyp.numbers <-
    match(hyp.nominal,
          c('Oral', 'Oral/Nasal', 'Nasal', 'Zygomatic', 'Orbit', 'Base', 'Vault'))

color2D.matplot(allo.data $ le.Wanc[order(hyp.numbers), order(hyp.numbers)])

allo.data $ Wanc2plot <- allo.data $ le.Wanc[order(hyp.numbers), order(hyp.numbers)]

colnames(allo.data $ Wanc2plot) <- rownames(allo.data $ Wanc2plot) <- NULL

allo.data $ Wanc2plot <- adply(allo.data $ Wanc2plot, 1:2)

colnames(allo.data $ Wanc2plot) <- c('row', 'col', 'value')

allo.data$Wanc2plot$trait.row <-
    factor(rownames(allo.data $ le.Wanc) [order(hyp.numbers)] [allo.data $ Wanc2plot $ row], levels = rownames(allo.data $ le.Wanc) [order(hyp.numbers)])

allo.data$Wanc2plot$trait.col <-
    factor(colnames(allo.data $ le.Wanc) [order(hyp.numbers)] [allo.data $ Wanc2plot $ col], levels = colnames(allo.data $ le.Wanc) [order(hyp.numbers)])

allo.data$Wanc2plot$hyp.row <-
    factor(hyp.nominal[order(hyp.numbers)] [allo.data $ Wanc2plot $ row],
           levels = c('Oral', 'Oral/Nasal', 'Nasal', 'Zygomatic',
            'Orbit', 'Base', 'Vault'))

allo.data$Wanc2plot$hyp.col <-
    factor(hyp.nominal[order(hyp.numbers)] [allo.data $ Wanc2plot $ col],
           levels = c('Oral', 'Oral/Nasal', 'Nasal', 'Zygomatic',
            'Orbit', 'Base', 'Vault'))

system('mkdir allo')

ggsave(
    'allo/Wancestral.pdf',
    ggplot(allo.data $ Wanc2plot) +
    geom_tile(aes(x = trait.col, y = trait.row, fill = value)) +
    scale_fill_viridis('Value', option = 'B', direction = -1) +
    scale_x_discrete(labels = hyp.nominal[order(hyp.numbers)]) +
    scale_y_discrete(limits = rev(levels(allo.data $ Wanc2plot $ trait.row))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab('Hypothesis') +
    ylab('Traits'),
    width = 12, height = 12)




