&record0
 input_file     = '%(mm5file)s' /

&record1
 start_year     =  %(syear)s
 start_month    =  %(smonth)s 
 start_day      =  %(sday)s
 start_hour     =  %(shour)s
 end_year       =  %(eyear)s                          
 end_month      =  %(emonth)s                          
 end_day        =  %(eday)s                         
 end_hour       =  %(ehour)s                           
 interval       = 10800 /

&record2
 pressure_bu_no_sfc_Pa = 100000 , 97500, 95000, 92500, 90000,  87500,
                         85000, 82500, 80000, 77500, 75000, 72500, 70000,
                         65000, 60000, 55000, 50000, 45000, 40000, 35000, 
                         30000, 25000, 20000, 15000, 10000 /
                                  

&record3
 print_info            = .false. /
