# main.py
import sys
try:
 x = int(input('Choose model (FOK:1, SOK:2, GOK:3): '))
 if x == 1:
    print('You have chosen the model FOK\nPlease go to the csv file and select the TL curve to analyze!')
    import FOK
 elif x == 2:
    print('You have chosen the model SOK\nPlease go to the csv file and select the TL curve to analyze!')
    import SOK
 elif x == 3:
    print('You have chosen the model GOK\nPlease go to the csv file and select the TL curve to analyze!')
    import GOK
 else:
    print('You have not selected any model, please go to main function to run it again')
except:
 print("Finish!")