library(rsconnect)

rsconnect::setAccountInfo(name='plantbreeding', 
                          token='6B041505BB404821E433D65CBB835BE9',
                          secret='eYYUzbJU/vHuEfzn/NM9CKQ/WI77e5FG9bRduVO6')
deployApp(account = 'plantbreeding', appName = "PBR_32803_qtl")

