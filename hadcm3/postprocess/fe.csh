#!/bin/csh

#foreach param (prc prw psl rlds rlus rlut rsdt rsus rsut snd snm tauu tauv ts zmla)
#set param="ta"
foreach param (cli clw hus)
 foreach level (1000 975 950 925 900 875 850 825 800 775 750 725 700 650 600 550 500 450 400 350 300 250 200 150 100)
  foreach scen (C)
    /usr/local/python/bin/python narccap_post.py $scen F $param $level
  end
 end
end
