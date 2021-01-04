# ----------------------------------------------------------------------------------------------------
# OUTIL ENVOI DE COURRIEL
# External packages
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : 22 avril 2020
# Last version  : 22 avril 2020
# ----------------------------------------------------------------------------------------------------

email.conv <- function(adresse.email,
                       pwd.server,
                       host.server,
                       port.server,
                       sujet="Message automatique de R",
                       message.email,
                       image.a.inclure=NULL){
  
  # Entête HTML pour corps du courriel
  html.content <- c(
    paste('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">',
          '<html xmlns="http://www.w3.org/1999/xhtml" lang="en-GB">',
          '<head>',
          '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />',
          '<title>R code RESULTS</title>',
          '<meta name="viewport" content="width=0.9*device-width, initial-scale=1.0"/>',
          '<style type="text/css">',
          '  td {',
          '    text-align: right;',
          '  }',
          '@media only screen and (max-width: 480px){',
          '  .dm-user-inserted-image{',
          '    height:auto !important;',
          '    max-width:600px !important;',
          '    width: 100% !important;',
          '  }',
          '}',
          '@media screen and (max-width:480px) {',
          '  table {',
          '    max-width: 100%!important;',
          '    height:auto;',
          '    display:block;',
          '  }',
          '  td {',
          '    width: 100%!important;',
          '  }', 
          '}',
          '</style>',
          '</head>',
          '</html>',
          '<body style="margin: 0cm; padding: 0.5cm;">',
          '<table style="float:right" border="1" cellpadding="0" cellspacing="0" width="100%">',
          '<tbody>',
          sep="\n",collapse="\n")
  )
  
  # Ajout des tables au corps du HTML
  for (message.part in message.email){
    html.content <- html.content %>% 
      append(
        paste0(toHTML(as.data.frame(message.part),
                      class.handlers = list(double = function(x) round(x,3))),
               '\n', collapse="\n"
        )
      )
  }
  # Clôture des tableaux
  html.content <- html.content %>% append(
    paste('</tbody>',
          '</table>',
          '</body>',
          sep="\n")
  )
  
  # Initialisation de l'enveloppe
  courriel <- envelope()
  
  # Ajout des ADRESSES COURRIELS à l'objet de courriel
  courriel <- courriel %>%
    from(adresse.email) %>%
    to(adresse.email) 
  
  # Ajout du SUJET du courriel
  courriel <- courriel %>% subject(sujet)
  
  # Ajout du CORPS du courriel
  courriel <- courriel %>% html(
    content=HTML(html.content,collapse="\n"),
    disposition = "inline",
    charset = "utf-8",
    encoding = "quoted-printable")
  
  # Ajout de l'heure au bas du tableau
  courriel <- courriel %>% 
    text(as.character(Sys.time()))
  
  # Inclusion d'une image (Ajout ggp)
  if (!is.null(image.a.inclure)){
    courriel <- courriel %>% attachment(image.a.inclure)
  }
  
  # Connection au serveur de messagerie
  smtp <- server(host = host.server,
                 port = port.server,
                 username = adresse.email,
                 password = pwd.server)
  
  # Effectue l'envoi final
  smtp(courriel, verbose = F)
}
