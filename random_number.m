%  generate the random number satisfied discrete distribute law
function random_num=random_number(value_choice,probability)
  len=length(value_choice) ;     
  pr=rand(1);
  nump=0;
     for j=1:len
          nump=nump+probability(j);
          if nump>=pr
             random_num=value_choice(j);
             return
          end
     end        
end 

       

   