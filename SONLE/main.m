//
//  main.m
//  SONLE
//
//  Created by Андрей Рычков on 16.03.14.
//  Copyright (c) 2014 Andrey Rychkov. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "MKSONLESolution.h"
#import "MKMatrix.h"

int main(int argc, const char * argv[])
{

    @autoreleasepool
    {
        NSNumber *eps = [NSNumber numberWithDouble:10e-12];
        NSArray *appr = [NSArray arrayWithObjects:
                         @0.5,
                         @0.5,
                         @1.5,
                         @-1.0,
                         @-0.5,
                         @1.5,
                         @0.5,
                         @-0.5,
                         @1.5,
                         @-1.5,
                         nil];
        
        
        MKSONLESolution *SONLE = [[MKSONLESolution alloc] initWithEpsilon:eps
                                                            approximation:appr
                                                        numberOfVariables:10];
        NSArray *x = [SONLE solveWithModified:YES index:5];
        
        for (NSNumber *number in x)
        {
            printf("\n %.12f  ", [number doubleValue]);
        }
        
        printf("\n \n");
    }
    return 0;
}
