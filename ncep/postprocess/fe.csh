#!/bin/csh

foreach param (wa clw cli hus)
#foreach param (ta ua va)
#foreach param (pr prc prw psl rlds rlus rlut rsdt rsus rsut snd snm tauu tauv ts zmla) 
 foreach level (1000 975 950 925 900 875 850 825 800 775 750 725 700 650 600 550 500 450 400 350 300 250 200 150 100)
    /usr/local/python/bin/python narccap_post_ncep.py C $param $level
  end
end

