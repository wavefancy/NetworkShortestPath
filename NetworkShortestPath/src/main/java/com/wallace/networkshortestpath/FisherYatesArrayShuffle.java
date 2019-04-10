/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.wallace.networkshortestpath;

import java.util.concurrent.ThreadLocalRandom;
import java.util.Arrays;
import java.util.List;

/**
 * Shuffle array by Fisher-Yates Algorithm.
 * https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * http://www.programming-algorithms.net/article/43676/Fisher-Yates-shuffle
 * -- To shuffle an array a of n elements (indices 0..n-1):
 * for i from n−1 downto 1 do
 *    j ← random integer such that 0 ≤ j ≤ i
 *    exchange a[j] and a[i]
 * @author wallace
 */
public class FisherYatesArrayShuffle {
    
    /**
     * Shuffle Array in place.
     * @param inArr 
     */
    public static void Shuffle(String[] inArr){
        ThreadLocalRandom random = ThreadLocalRandom.current();
        for (int i = inArr.length -1 ; i > 0; i--) {
            int index = random.nextInt(i+1);
            String t = inArr[i];
            inArr[i] = inArr[index];
            inArr[index] = t;
        }
    }
    
    /**
     * Random pick length elements from an input array(inArr), length <= inArr.length - 1.
     * Thread saft, only one thread can access this code at a particuar time point.
     * @param inArr
     * @param length
     * @return 
     */
    public static synchronized List<String> Shuffle(String[] inArr, int length){
        ThreadLocalRandom random = ThreadLocalRandom.current();
        if (length >= inArr.length) {
            System.err.println("ERROR IN: shortestPath.edu.princeton.cs.algs4.FisherYatesArrayShuffle.Shuffle()");
            System.err.println("new list should be as less or equal length as input array.");
            System.exit(-1);
        }
        for (int i = 0; i < inArr.length; i++) {
            int index = random.nextInt(i,inArr.length);
            String t = inArr[i];
            inArr[i] = inArr[index];
            inArr[index] = t;
            if (i == length) {
                break;
            }
        }
        
        //return  Arrays.asList(inArr).subList(0, length).toArray();
        return Arrays.asList(Arrays.copyOf(inArr, length));
    }
   

    public static void main(String[] args) {
        String[] temp = {"1","2","3","4"};
        //FisherYatesArrayShuffle.Shuffle(temp);
        for (String string : FisherYatesArrayShuffle.Shuffle(temp, 4)) {
            System.out.println(string);
        }
    }
}
