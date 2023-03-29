;;;;;;;;;; FIRST INDUSTRIAL REVOLUTION - NEW VERSION (?!) - 17/01/21

;set random-seed 100

globals [
  ;used to monitor the single runs
  aggregate-farm-production
  aggregate-industry-production
  aggregate-services
  aggregate-metabolism
  aggregate-goods-demand
  aggregate-services-demand
  aggregate-profits
  aggregate-dividends
  average-capital-firms
  average-capital-farms
  average-profits-firms
  average-profits-farms
  mean-profits-farms
  mean-profit-firms
  yearly-interest

  ;counting
  industrialized-firms
  failed-firms
  new-firms
  local-number
  unemployed
  employed
  bonds-sold
  kids-number
  number-of-firms2
  number-of-workers2
  initial-food-prices2
  initial-good-prices2
  initial-service-prices2


  ;means
  mean-farm-price
  mean-firm-price
  mean-service-price
  mean-salaries
  mean-salaries-firm-workers
  mean-salaries-farm-workers
  mean-salaries-public-workers

  ;checks for avoding bugs
  industrial-switching-check
  labor-check
  labor-check-1
  labor-check-2
  production-check
  demand-check
  farms-market-check
  firms-market-check
  service-market-check
  dividends-check
  price-adjusting-check
  agent-replacement-check
  government-check-1
  government-check-2
  government-check-3
  kill-switch

  ;actual macro variables
  taxes
  GDP-income
  GDP-spending
  real-GDP-spending
  goods-income-value
  labor-income-value
  land-income-value
  service-income-value
  profit-income-value
  inflation
  exports
  imports

  ;initial parameters, transferred to choosers
;  number-of-cities
;  number-of-coal
;  initial-labor-price-workers
;  initial-labor-price-employers
;  initial-household-wealth
;  bourgeoisie-ratio-wealth
;  nobles-ratio-wealth
;  initial-industrial-productivity
;  firm-productivity
;  farm-productivity
;  firm-labor ;max
;  firm-industrial-labor ;max when switching
;  service-productivity
;  initial-food-price
;  initial-goods-price
;  initial-service-price
;  initial-labor-price
;  initial-capital
;  initial-capital
;  substitute-capital
;  government-initial-wealth
;  number-of-workers
;  number-of-bourgeoisie
;  number-of-nobles
;  number-of-firms
;  number-of-farms

  ;; other general parameters for setup
;  distribution ;decides range of distribution
;  reproduction ;how high the value of kids has to be to spawn new agent
;  markup ;how much markup compared to productivity
;  distance-setup ;how far can they see
;  tax-rate ; the ratio of taxation, usally very low
;  price-change-chance ;probability of price change, the higher it is the more likely it will be
;  price-delta ;how much can the price change

  ;industrial transition
;  trans-cost ;how much should it cost
;  safe-zone ;how much more capital than cost for industrial switching
;  industrial-switch-probability ;the higher it is the less likely it will happen
;  time ;when industrial revolution

;  ;labor-market
;  percent-government-workers ;general percentage of the labor force that the government should hire

  ;new substitution mechanism?
  failed-tick
  starting-seed 
]


;;;; BREEDS
breed [workers worker]
breed [bourgeoisie bourgeois]
breed [nobles noble]
breed [firms firm]
breed [farms farm]
breed [government gov]
breed [foreign-market foreign]

undirected-link-breed [ownership own]


;;; BREEDS VARIABLES
workers-own [
  previous-wealth
  wealth ;available wealth
  metabolism ;demand for food
  needs ; demand for goods
  desires; demand for service
  kids; if satisfied the agent will reproduce
  unemployed?; if employed or not
  salary ;demanded salary
  firm? ;hired by a firm
  farm? ;hired by a farm
  state? ;hired by the state
  ]

bourgeoisie-own [
  previous-wealth
  wealth ;available wealth
  metabolism ;demand for food
  needs ; demand for goods
  desires; demand for service
  kids; if satisfied the agent will reproduce
  bonds  ;number of owned bonds
  owner? ;if they own a firm
  price ;initial service price
  productivity ;service productivity
  service; services produced
]

nobles-own [
  previous-wealth
  wealth ;available wealth
  metabolism ;demand for food
  needs ; demand for goods
  desires; demand for service
  bonds ;number of owned bonds
  kids
]

firms-own [
  capital  ;available funds
  previous-capital ;needed for profit calculation
  stock ;unsold past products
  labor ;workers hired this turn
  coal? ;binary variable for the presence of nearby coal
  failed?
  industrywannabe?
  industry?
  max-labor
  max-industrial-labor
  productivity
  industrial-productivity
  price
  labor-price
  profits
  failed?
]

farms-own[
  capital
  previous-capital
  stock
  labor
  failed?
  max-labor
  productivity
  price
  labor-price
  profits
  failed?
]

government-own[
  wealth
  labor
  labor-price
  metabolism
  needs
  desires
  interest-rate
  wage
  bonds
]

foreign-market-own[
  wealth
  metabolism
  needs
  desires
  goods
  food
  food-price
  goods-price
]

patches-own [
  city?
  coal2?
]

to setup
  clear-all
  clear-globals
  set starting-seed new-seed
  random-seed starting-seed
 ;;; parameters initial setup, remember to change for experiment
  setup-patches
  setup-workers
  setup-bourgeoisie
  setup-nobles
  setup-firms
  setup-farms
  if government-features = true [setup-government]
  if government-features = false [set tax-rate 1]
  if foreign-market-features = true [setup-foreign-market]
  set initial-food-prices2 mean [price] of farms
  set initial-good-prices2 mean [price] of firms
  set initial-service-prices2 mean [price] of bourgeoisie
  reset-ticks
end

to setup-patches
  ask patches [
    set pcolor green
    set city? false
    set coal2? false
  ]
  let fields patches with [pcolor = green]
  ask n-of number-of-cities fields [
    set city? true
    set pcolor yellow
  ]
    set fields patches with [pcolor = green]
    ask n-of number-of-coal fields [
      set coal2? true
      set pcolor black
    ]
end

to setup-workers
  create-workers number-of-workers [
    set color red
    set size 0.1
    set shape "person"
    setxy random-xcor random-ycor
    set wealth random-normal initial-household-wealth (initial-household-wealth / distribution)
    set previous-wealth wealth
    set salary random-normal initial-labor-price (initial-labor-price / distribution)
    set kids random-normal (reproduction / 2) (reproduction / distribution)
  ]
end

to setup-bourgeoisie
  create-bourgeoisie number-of-bourgeoisie [
    set color orange
    set size 0.1
    set shape "person"
    setxy random-xcor random-ycor
    set wealth random-normal (initial-household-wealth * bourgeoisie-ratio-wealth) (initial-household-wealth * bourgeoisie-ratio-wealth / distribution)
    set price random-normal initial-service-price (initial-service-price / distribution)
    set kids random-normal (reproduction / 2) (reproduction / distribution)
    set productivity round random-normal service-productivity (service-productivity / distribution)
    set previous-wealth (wealth / 5)
    set owner? false
  ]
  let cities patches with [pcolor = yellow]
  let bourgeoisie-per-city ( round (((count bourgeoisie ) / 2)  / (count cities)))
  let bourgeoisie-moving n-of ((count bourgeoisie) / 2) bourgeoisie
  ask cities [
    ask n-of bourgeoisie-per-city bourgeoisie-moving [
      move-to myself
      right random 360
      forward random-float 1
    ]
  ]


end

to setup-nobles
  create-nobles number-of-nobles [
    set color blue
    set size 0.1
    set shape "person"
    setxy random-xcor random-ycor
    set wealth random-normal (initial-household-wealth * nobles-ratio-wealth) (initial-household-wealth * nobles-ratio-wealth / distribution)
    set previous-wealth (wealth / 10)
  ]
end

to setup-firms
  let owners n-of number-of-firms bourgeoisie
  ask owners [
    set owner? true
    hatch-firms 1 [
      create-own-with myself
      set shape "house"
      set color orange
      set size 0.2
      right random 360
      forward random-float 0.5
      set coal? false
      set industrywannabe? false
      set industry? false
      set failed? false
      set profits 0
      set capital random-normal initial-capital (initial-capital / distribution)
      set productivity random-normal firm-productivity (firm-productivity / distribution)
      set max-labor random-normal firm-labor (firm-labor / distribution)
      set max-industrial-labor random-normal firm-industrial-labor (firm-industrial-labor / distribution)
      set industrial-productivity random-normal initial-industrial-productivity (initial-industrial-productivity / distribution)
      set labor-price random-normal initial-labor-price (initial-labor-price / distribution)
      set price ((labor-price / productivity) + random-float markup )
      let close in-radius distance-setup patches
      if any? close with [pcolor = black][set coal? true]
    ]
  ]
  set number-of-firms2 number-of-firms
end

to setup-farms
  let total-needs ((count workers) + ((count bourgeoisie) * 5) + ((count nobles) * 20))
  let pastures patches with [pcolor = green]
  let total-pastures count pastures
  ask pastures [
    sprout-farms 1 [
      let owner min-one-of nobles [distance myself]
      ask owner [create-own-with myself]
      set shape "house"
      set color violet
      set size 0.2
      set failed? false
      set profits 0
      set capital random-normal initial-capital (initial-capital / distribution)
      set labor-price random-normal initial-labor-price (initial-labor-price / distribution)
      set productivity 3
      set max-labor (round ((total-needs / total-pastures) / 3))
      set price ((labor-price / productivity) + random-float markup)
    ]
  ]
end

  to setup-government
  create-government 1 [
    set size 0.2
    set color black
    set interest-rate 0.035
    set shape "person"
    set wage initial-labor-price
    set wealth government-initial-wealth
    set metabolism 0
    set desires 0
    set wealth 0
  ]
end


  to setup-foreign-market
  let mean-food-price mean [price] of farms
  let mean-goods-price mean [price] of firms

  create-foreign-market 1
  ask foreign-market [
    set shape "person"
    set size 0.1
    set color yellow
    set metabolism 0
    set desires 0
    set wealth 0
    set food 0
    set food-price 0
    set goods 0
    set goods-price 0
  ]
end


to go
  set taxes 0
  set GDP-income 0
  set GDP-spending 0
  set real-GDP-spending 0
  set goods-income-value 0
  set labor-income-value 0
  set land-income-value 0
  set service-income-value 0
  set profit-income-value 0


;  set inflation 0
;  set exports 0
;  set imports 0
  if ticks > time [run-industrial-switching]
  run-labor-matching
  run-production
  run-demand
  run-food-market
  run-goods-market
  run-service-market
  run-dividends
  run-price-adjustment
  if government-features = true [run-government]
  run-substitution
  tick
end


to run-industrial-switching
  let mean-labor-price mean [labor-price] of firms
  ;let trans-cost2 (mean-labor-price * trans-cost-multiplier)
  set kill-switch 0
  set industrial-switching-check 0
  let not-industrialized firms with [industry? = false]
  let industrialized firms with [industry? = true]
  let coal-industry not-industrialized with [coal? = true]
  let wannabe-industries coal-industry with [industrywannabe? = true]
  set industrial-switching-check (industrial-switching-check + 1)
  let different (mean-labor-price - initial-labor-price)
  let cost-increase (different / initial-labor-price)
  ask wannabe-industries [
    set industrial-switching-check (industrial-switching-check + 1)
    set kill-switch (kill-switch + 1)
    if kill-switch > (number-of-firms * 2) [stop]
    let available-capital 0
    let safe-cost 0
    let number random-float 1
    let owner-capital 0
    let close-firms industrialized in-radius (distance-setup / 2)
    let number-industry min (list 5 ((count close-firms) / 2))
    set safe-cost ((1 + (random-float safe-zone)) * (trans-cost * cost-increase * (1 - (number-industry) / 10)))
    ;set safe-cost ((1 + (random-float safe-zone)) * (trans-cost2 * (1 - (number-industry) / 10)))
    set owner-capital [wealth] of one-of out-own-neighbors
    set available-capital (owner-capital + capital)
    if available-capital >= safe-cost and (number < industrial-switch-probability) [
      set industry? true
      set industrywannabe? false
      set industrialized-firms (industrialized-firms + 1)
      set max-labor max-industrial-labor
      set productivity industrial-productivity
      let early-capital capital
      set capital (capital + owner-capital - safe-cost)
      set capital (capital / 2)
      let payed capital
      ask out-own-neighbors [ set wealth payed ] ;owner and firm share the cost
    ]

  ]
end


;new model is unable to sustain revolution with no additional demand

to run-labor-matching
  set kill-switch 0
  set labor-check-1 0
  set labor-check-2 0
  set labor-check false
  ask workers [
    set unemployed? true
    ;set employed? false
    set state? false
    set firm? false
    set farm? false
    set labor-check true
    set previous-wealth wealth
  ]
  ;this exist to avoid problem when you have more workers than initial number
  set number-of-workers2 count workers
  set unemployed 0
  set employed 0
  set GDP-income 0
  set labor-income-value 0
  ; set labor-check 0
  ;Public Employment
  if government-features = true [
    let government-workers round random-normal (percent-government-workers * number-of-workers2) ( (percent-government-workers * number-of-workers2) / distribution)
    if government-workers < 0 [set government-workers (number-of-workers2 / 10)] ;safeguard
    let public-workforce n-of government-workers workers
    let public-salary [wage] of one-of government ;of one-of employed to prevent bug where wage is treated not as a value but as a list, resulting in error
    ask public-workforce [
      set wealth (wealth + public-salary)
      set unemployed? false
      set state? true
      set employed (employed + 1)
      set GDP-income (GDP-income + public-salary)
      set labor-income-value ( labor-income-value + public-salary)
      set labor-check-1 (labor-check-1 + 1)
    ]
    ask government [
      let labor-expense (wage * government-workers)
      set wealth (wealth - labor-expense)
      set labor count workers with [state? = true]
    ]
  ]

  ;;other employers
  let employers turtles with [shape = "house"]
  let actually-unemployed workers with [unemployed? = true]
  let average-food-production mean [max-labor] of farms
  let seasonality round random-normal 0 (average-food-production)
  ask employers [
    set previous-capital capital
    set labor 0
  ]
  let actual-employers employers with [ capital > labor-price]
  ifelse (any? actual-employers) and (any? actually-unemployed) [
    ask actual-employers [
      set kill-switch (kill-switch + 1)
      if kill-switch > (number-of-workers) [stop]
      let hirable round (capital / labor-price)
      ifelse color = violet [
        set hirable min list (max-labor +(random seasonality)) hirable
      ]
      [set hirable min list max-labor hirable]
      let offered-wage labor-price
      let local-unemployed actually-unemployed in-radius distance-setup with [salary <= offered-wage]
      if ((count local-unemployed) < 0) or (hirable < 0) [stop]
      if any? local-unemployed  [
        let hired up-to-n-of hirable local-unemployed
        ;up-to-n-of allows to pick less than demanded if the subsample is not big enough
        set labor count hired
        set capital (capital - (labor-price * labor))
        if any? hired [
          ask hired[

            if [color] of myself = orange [set firm? true]
            if [color] of myself = violet [set farm? true]
            set unemployed? false
            set wealth (wealth + offered-wage)
            set employed (employed + 1)
            set GDP-income (GDP-income + offered-wage)
            set labor-income-value ( labor-income-value + offered-wage)
            set actually-unemployed actually-unemployed with [self != myself]
            set labor-check-2 (labor-check-2 + 1)
          ]
        ]
      ]
        set actual-employers employers with [self != myself]
        stop
    ]
  ][stop]
  let unemp workers with [unemployed? = true]
  set unemployed count unemp
  let employed2 workers with [unemployed? = false]
  let public-employed employed2 with [state? = true]
  if any? employed2 [
    set mean-salaries (mean [salary] of employed2)
    let farm-workers employed2 with [farm? = true]
    let firms-workers employed2 with [firm? = true]
    if any? firms-workers [
      set mean-salaries-firm-workers (mean [salary] of employed2 )
    ]
    if any? farm-workers [
      set mean-salaries-farm-workers (mean [salary] of  farm-workers)
    ]
    if ((government-features = true) and (any? public-employed)) [
      set mean-salaries-public-workers (mean [salary] of public-employed)
    ]
  ]
end

to run-production
  set kill-switch 0
  set aggregate-farm-production 0
  set aggregate-industry-production 0
  set aggregate-services 0
  set production-check 0
  let producers turtles with [shape = "house"]
  ask producers [
    set kill-switch (kill-switch + 1)
    if kill-switch > (number-of-workers) [stop]
    set production-check (production-check + 1)
    set stock round ((labor * productivity) + stock)
    ifelse [color] of self = violet [
      set aggregate-farm-production (aggregate-farm-production + (labor * productivity)) ] [
     set aggregate-industry-production (aggregate-industry-production + (labor * productivity)) ]
    ]
  let service-bourgeoisie bourgeoisie with [owner? = false]
   ask service-bourgeoisie [
    ;if they have a firm they should not sell services
      set service round ( random-normal productivity (productivity / distribution))
      set aggregate-services ( aggregate-services + service)
  ]
end

to run-demand
  set kill-switch 0
  set mean-farm-price mean [price] of farms
  set mean-firm-price mean [price] of firms
  set mean-service-price mean [price] of bourgeoisie

  ask workers [
    set kill-switch (kill-switch + 1)
    if kill-switch > (number-of-workers * 2) [stop]
    let income (wealth - previous-wealth)
    let spending random-normal income (income / distribution)
    ifelse income > 0 [
      let food-propensity round (( spending / 2) / mean-farm-price )
      let goods-propensity round ((spending / 2) / mean-firm-price)
      if (food-propensity > 1) and (goods-propensity > 1) [
        set metabolism max list 1 food-propensity
        set needs max list 1 goods-propensity
      ]
    ] [
      set metabolism 1
      set needs 1
    ]
  ]

  set kill-switch 0

  ask bourgeoisie [
    set kill-switch (kill-switch + 1)
    if kill-switch > (number-of-bourgeoisie * 3) [stop]
    let income (wealth - previous-wealth)
    ifelse income > 0 [
      let spending random-normal income (income / distribution)
      let food-propensity round ((spending / 5) / mean-farm-price)
      let goods-propensity round ((spending / 4) / mean-firm-price)
      let service-propensity round ((spending / 4) / mean-firm-price)
      if (food-propensity > 1) and (goods-propensity > 1) [
        set metabolism max list 1 food-propensity
        set needs max list 1 goods-propensity
        set desires max list 1 service-propensity
      ]
    ][
      set metabolism 1
      set needs 1
    ]
  ]
  set kill-switch 0

  ask nobles [
    set kill-switch (kill-switch + 1)
    if kill-switch > (number-of-nobles *  3) [stop]
    let income (wealth - previous-wealth)
    ifelse income > 0 [
      let spending random-normal income (income / distribution)
      let food-propensity round ((spending / 8) / mean-farm-price)
      let goods-propensity round ((spending / 8 * 3 ) / mean-firm-price)
      let service-propensity round ((spending / 8 * 3) / mean-firm-price)
      if (food-propensity > 1) and (goods-propensity > 1) [
        set metabolism max list 1 food-propensity
        set needs max list 1 goods-propensity
        set desires max list 1 service-propensity
      ]
    ][
      set metabolism 1
      set needs 1
    ]
  ]

  if government-features = true [
    ask government [
      ;set metabolism random-normal (aggregate-farm-production * government-demand / 2) ((aggregate-farm-production * government-demand / 2) / distribution)
      set needs random-normal  round (aggregate-industry-production * government-demand / 2) ((aggregate-industry-production * government-demand / 2) / distribution)
      set desires random-normal round (aggregate-services * government-demand ) ((aggregate-services * government-demand ) / distribution)
    ]
  ]




end



to run-food-market
  set kill-switch 0
  let consumers turtles with [ shape = "person"]
  set aggregate-metabolism sum [metabolism] of consumers
  let available-food aggregate-farm-production
  set land-income-value 0
  set goods-income-value 0
  set service-income-value 0
  set farms-market-check 0
  set firms-market-check 0
  set service-market-check 0
  set government-check-1 0
  set government-check-2 0
  set government-check-3 0



  ;; food market
  let mean-price mean [price] of farms
  let minimum-price min [price] of farms
  let demand consumers with [ (wealth > mean-price) and (metabolism > 0)]
  let supply farms with [stock > 0]

  while [any? demand] [
    ask demand [
      let close-supply 0
      ifelse [color] of self = black [
          set close-supply supply in-radius distance-setup ][
        set close-supply supply
        set government-check-1 (government-check-1 + 1)
      ]
      while [(metabolism > 0) and (wealth > minimum-price) and ( any? close-supply)] [
        set kill-switch (kill-switch + 1)
        if kill-switch > (number-of-workers * 5) [stop]
        ;let black-demand demand with [color = black]
        ;if any? black-demand [ask black-demand [set close-supply supply]]
        let low-price min-one-of close-supply [price]
        ;this allows for the price of low-price to be seen by the consumer
        let chosen-price [price] of low-price
        ;this choses the lowest value between the available funds of the consumer and its desire to consume
        let max-buying min list (round (wealth / (chosen-price * tax-rate))) (metabolism)
        ;this allows the consumer to see the stock of the supplier
        let available-stock [stock] of low-price
        ;this choses the lowest value between the desired quantity by the consumer and the available quantity at the seller
        let buying min list max-buying available-stock
        let spending (buying * (chosen-price * tax-rate))
        set wealth (wealth - spending )
        set land-income-value (land-income-value + spending)
        set GDP-spending (GDP-spending + spending)
        set taxes (taxes + (spending - (chosen-price * buying)))
        set farms-market-check (firms-market-check + 1)
          ask low-price [
          set capital (capital + spending)
          set stock (stock - buying)
          set real-GDP-spending (real-GDP-spending + (buying * initial-food-prices2))

          ;if stock is depleted seller is removed from supply pool
          if stock = 0 [
            set supply supply with [self != myself]
            set close-supply close-supply with [self != myself]
            set farms-market-check (farms-market-check + 1)
          ]
        ]
        ;  set demand consumers with [ (wealth > mean-price) and (needs > 0)]
        ifelse buying = metabolism [
          set kids (kids + 1)
          set metabolism 0
          set demand demand with [self != myself]
        ] [set metabolism (metabolism - buying)
        ]
        let lower-price [price] of low-price
        if (metabolism = 0) or ( wealth < lower-price) [
          set demand demand with [self != myself]
          stop
        ]
      ]
    ]
    stop
  ]
end

 to run-goods-market
  set kill-switch 0
  let consumers turtles with [ shape = "person"]
  set aggregate-goods-demand sum [needs] of consumers
  let available-goods aggregate-industry-production
  let mean-price mean [price] of firms
  let minimum-price min [price] of firms
  let demand consumers with [ (wealth > mean-price) and (needs > 0)]
  let supply firms with [stock > 0]
  ;; GOODS MARKET
  while [any? demand] [
    ask demand [
      let close-supply 0
      ifelse [color] of self = black [
        set close-supply supply in-radius distance-setup ][
        set close-supply supply
        set government-check-1 (government-check-1 + 1)
      ]
      while [(needs > 0) and (wealth > minimum-price) and ( any? close-supply)] [
        set kill-switch (kill-switch + 1)
        if kill-switch > (number-of-workers * 5) [stop]
        ;let black-demand demand with [color = black]
        ;if any? black-demand [ask black-demand [set close-supply supply]]
        let low-price min-one-of close-supply [price]
        ;this allows for the price of low-price to be seen by the consumer
        let chosen-price [price] of low-price
        ;this choses the lowest value between the available funds of the consumer and its desire to consume
        let max-buying min list (round (wealth / (chosen-price * tax-rate))) (needs)
        ;this allows the consumer to see the stock of the supplier
        let available-stock [stock] of low-price
        ;this choses the lowest value between the desired quantity by the consumer and the available quantity at the seller
        let buying min list max-buying available-stock
        let spending (buying * (chosen-price * tax-rate))
        set wealth (wealth - spending )
        set goods-income-value (goods-income-value + spending)
        set GDP-spending (GDP-spending + spending)
        set taxes (taxes + (spending - (chosen-price * buying)))
        set firms-market-check (firms-market-check + 1)
        set government-check-2 ( government-check-2 + 1)
        ask low-price [
          set capital (capital + spending)
          set stock (stock - buying)
          set real-GDP-spending (real-GDP-spending + (buying * initial-good-prices2))
          ;if stock is depleted seller is removed from supply pool
          if stock = 0 [
              set supply supply with [self != myself]
            set close-supply close-supply with [self != myself]
            set firms-market-check (firms-market-check + 1)
            ]
        ]
        ;  set demand consumers with [ (wealth > mean-price) and (needs > 0)]
        let lower-price [price] of low-price
        if (needs = 0) or ( wealth < lower-price) [
          set demand demand with [self != myself]
          stop
        ]
      ]
    ]
    stop
  ]
end

to run-service-market
  set kill-switch 0
  let consumers turtles with [ shape = "person"]
  set aggregate-services-demand sum [desires] of consumers
  let available-services aggregate-services
  let mean-price mean [price] of bourgeoisie
  let minimum-price min [price] of bourgeoisie
  let demand consumers with [ (wealth > mean-price) and (desires > 0) ]
  let supply bourgeoisie with [(service > 0)] ; and (owner? = false)
                                              ;; SERVICE MARKET
  while [any? demand] [
    ask demand [
      let close-supply 0
      ifelse [color] of self = black [
        set close-supply supply in-radius (distance-setup * service-distance-multiplier) ][
        set close-supply supply
        set government-check-1 (government-check-1 + 1)
      ]
      while [(desires > 0) and (wealth > minimum-price) and ( any? close-supply)] [
        set kill-switch (kill-switch + 1)
        if kill-switch > (number-of-workers * 5) [stop]
        ;let black-demand demand with [color = black]
        ;if any? black-demand [ask black-demand [set close-supply supply]]
        let low-price min-one-of close-supply [price]
        ;this allows for the price of low-price to be seen by the consumer
        let chosen-price [price] of low-price
        ;this choses the lowest value between the available funds of the consumer and its desire to consume
        let max-buying min list (round (wealth / (chosen-price * tax-rate))) (desires)
        ;this allows the consumer to see the stock of the supplier
        let available-stock [service] of low-price
        ;this choses the lowest value between the desired quantity by the consumer and the available quantity at the seller
        let buying min list max-buying available-stock
        let spending (buying * (chosen-price * tax-rate))
        set wealth (wealth - spending )
        set service-income-value (service-income-value + spending)
        set GDP-spending (GDP-spending + spending)
        set taxes (taxes + (spending - (chosen-price * buying)))
        ask low-price [
          set wealth (wealth + spending)
          set service (service - buying)
          set real-GDP-spending (real-GDP-spending + (buying * initial-service-prices2))
          set service-market-check (service-market-check + 1)
          ;if stock is depleted seller is removed from supply pool
          if service = 0 [
            set supply supply with [self != myself]
            set close-supply close-supply with [self != myself]
          ]
        ]
        ;  set demand consumers with [ (wealth > mean-price) and (needs > 0)]
        let lower-price [price] of low-price
        if (needs = 0) or ( wealth < lower-price) [
          set demand demand with [self != myself]
          stop
        ]
      ]
    ]
    stop
  ]
end




to run-dividends
  set kill-switch 0
  set aggregate-dividends 0
  set aggregate-profits 0
  set failed-tick 0

    let households turtles with [(color = orange) or (color = "blue") and (shape = "person") ]
  ask households [
    set previous-wealth wealth
  ]

  ;; farms
  ask farms [
    if capital <= mean-salaries [
      ask out-own-neighbors [
        let transfer min ( list average-capital-farms ( 0.5 * wealth))
        set wealth (wealth - transfer)
        ask myself [set capital (capital + transfer)]
      ]
    ]
  ]


  ;; firms
    ask firms [
    if capital < labor-price [
      set failed-firms ( failed-firms + 1)
      set failed-tick (failed-tick + 1)
      ask out-own-neighbors [
        let owned count out-own-neighbors
        if owned < 1 [set owner? false]
      ]
      die
    ]
  ]


  let employers turtles with [shape = "house"]
  ask employers [
    set profits (capital - previous-capital)
  ]
  let profitable employers with [profits > 0]
  if any? profitable [
    ask profitable [
      set aggregate-profits (aggregate-profits + profits)
    ; set GDP-income (GDP-income + profits) ?
      let min-dividends (profits / 2)
      let random-dividends (random (profits / 2))
      let dividends (min-dividends + random-dividends)
      if color = orange [
        if capital > (trans-cost * (safe-zone * 2)) [
          set dividends ( capital - (trans-cost * (safe-zone * 2)))
        ]
      ]
      if color = violet [
        let limit (max-labor * labor-price)
        if capital > limit [ set dividends (dividends + (capital - limit))]
        ]
      set aggregate-dividends (aggregate-dividends + dividends)
      set capital (capital - dividends)
      let owners count out-own-neighbors ;;just in case a family splits and both mantain ownership
      ask out-own-neighbors [
        set wealth (wealth + (dividends / owners))
        set profit-income-value (profit-income-value + dividends)
      ]
    ]
  ]
  set average-capital-firms mean [capital] of firms
  set average-capital-farms mean [capital] of farms
  set mean-profits-farms mean [profits] of farms
  set mean-profit-firms mean [profits] of firms
  ask firms [
    if stock = 0 [set industrywannabe? true]
  ]
end

to  run-price-adjustment
  let employers turtles with [shape = "house"]
  ask employers [
    let labor-number (previous-capital / labor-price)
    let random-number random-float 1
    if (labor < labor-number) and (random-number > price-change-chance) [
      set labor-price ( labor-price + (random-normal price-delta (price-delta / distribution)))
    ]
    if (labor = labor-number) and (random-number > price-change-chance) [
      set labor-price ( labor-price - (random-normal price-delta (price-delta / distribution)))
      ]
  ]


  ask workers [
    let mean-food-price mean [price] of farms
    let random-number random-float 1
    if (unemployed? = false) and (random-number < price-change-chance) [
      set salary ( salary + (random-normal price-delta (price-delta / distribution))) ]
    if (unemployed? = true) and (random-number < price-change-chance) [
      set salary ( salary - (random-normal price-delta (price-delta / distribution)))]
    if salary < mean-food-price [set salary (mean-food-price + random-float 1)]
  ]



  ;market prices
  ask employers [
    let random-number random-float 1
    if (stock = 0) and (random-number < price-change-chance) and (labor > 0) [
     set price ( price + (random-normal price-delta (price-delta / distribution)))
    ]
    if (stock > 0) and (random-number < price-change-chance) [
      set price ( price - (random-normal price-delta (price-delta / distribution)))
    ]
    if (labor-price / productivity) > price  [
      set price ((labor-price / productivity) + random-float markup)
    ]
  ]

  ask bourgeoisie [
    let random-number random-float 1
    if (service = 0) and (random-number < price-change-chance) [
      set price ( price + (random-normal price-delta (price-delta / distribution)))
    ]
    if (service > 0) and (random-number < price-change-chance) [
      set price ( price - (random-normal price-delta (price-delta / distribution)))
    ]
    set price-adjusting-check (price-adjusting-check + 1)
    if price < 1 [set price 1]
  ]

  set mean-farm-price mean [price] of farms
  if any? firms [ set mean-firm-price mean [price] of firms]
  set mean-service-price mean [price] of bourgeoisie

end

to  run-government
  set kill-switch 0
  ask government [
    set bonds-sold 0
    set bonds 0
    set wealth (wealth + taxes)
    set taxes 0
    ;;the following split is necessary as not all agents have the bond variable
    let households turtles with [shape = "person"]
    let possible-bond-holders households with [(color = blue) or (color = orange)]
    let bond-holders possible-bond-holders with [bonds > 0]
    let effective-interests (1 + interest-rate)
    let interest-paid ((sum [bonds] of bond-holders) * effective-interests)
    ask bond-holders [
      set wealth (wealth + (bonds * effective-interests))
    ]
    set wealth (wealth - interest-paid)
    set yearly-interest interest-paid
    if wealth < 0  [
      let value abs wealth
      set bonds round ((value / 100) - 1)
      let possible-bond-buyers households with [(color = blue) or (color = orange)]
      let bond-buyers possible-bond-buyers with [wealth > 100]
      while [(bonds > 0) and (any? bond-buyers)] [
        set kill-switch (kill-switch + 1)
        if kill-switch > (number-of-workers * 5) [stop]
        ask bond-buyers [
          set bonds (bonds + 1)
          set wealth (wealth - 100)
          ask myself [
            set bonds (bonds - 1)
            set wealth (wealth + 100)]
          if wealth < 100 [
            set bond-buyers bond-buyers with [self != myself]
          ]
        ]
      ]
    ]
  ]

end

to run-substitution
  set kill-switch 0
  let possible-breeders turtles with [(color = red) or (color = orange) and (shape = "person") ]
  let breeders possible-breeders with [  (kids > reproduction)]
  ask breeders  [
    hatch 1 [
      set wealth (wealth / 2)
      set kids 0
      ]
    set wealth (wealth / 2)
    set kids-number (kids-number + 1)
    set kids 0
  ]

  set average-capital-firms mean [capital] of firms
  let substitute-capital random-normal average-capital-firms (average-capital-firms / distribution)
  let mean-labor-price mean [labor-price] of firms
  let possible-replacement bourgeoisie with [wealth > substitute-capital]
  let existing-firms count firms
  let mean-productivity mean [productivity] of firms
  let mean-industrial-productivity mean [industrial-productivity] of firms
  let missing-firms (number-of-firms2 - (count firms))
  if missing-firms >= 0 [
    ask up-to-n-of missing-firms possible-replacement [
      hatch-firms 1 [
        create-own-with myself
        set shape "house"
        set color orange
        set size 0.1
        right random 360
        forward random-float 0.5
        set capital substitute-capital
        set coal? false
        set industrywannabe? false
        set industry? false
        set failed? false
        set profits 0
        set capital random-normal initial-capital (initial-capital / distribution)
        set productivity random-normal mean-productivity (mean-productivity / distribution)
        set max-labor random-normal firm-labor (firm-labor / distribution)
        set max-industrial-labor random-normal firm-industrial-labor (firm-industrial-labor / distribution)
        set industrial-productivity random-normal mean-industrial-productivity (mean-industrial-productivity / distribution)
        set labor-price random-normal  mean-labor-price ( mean-labor-price / distribution)
        set price ((labor-price / productivity) + random-float markup )
        set new-firms (new-firms + 1)
        let close in-radius distance-setup patches
        if any? close with [pcolor = black][set coal? true]
        ask myself [
          set wealth (wealth - substitute-capital)
          set owner? true
        ]
      ]
    ]
  ]
  if missing-firms = 0[
    let more-people round (kids-number / 10)
    let difference ( more-people + number-of-firms2 - existing-firms)
      ask up-to-n-of difference possible-replacement [
        hatch-firms 1 [
create-own-with myself
        set shape "house"
        set color orange
        set size 0.1
        right random 360
        forward random-float 0.5
        set capital substitute-capital
        set coal? false
        set industrywannabe? false
        set industry? false
        set failed? false
        set profits 0
        set capital random-normal initial-capital (initial-capital / distribution)
        set productivity random-normal mean-productivity (mean-productivity / distribution)
        set max-labor random-normal firm-labor (firm-labor / distribution)
        set max-industrial-labor random-normal firm-industrial-labor (firm-industrial-labor / distribution)
        set industrial-productivity random-normal mean-industrial-productivity (mean-industrial-productivity / distribution)
        set labor-price random-normal  mean-labor-price ( mean-labor-price / distribution)
        set price ((labor-price / productivity) + random-float markup )
        set new-firms (new-firms + 1)
        let close in-radius distance-setup patches
        if any? close with [pcolor = black][set coal? true]
        ask myself [
          set wealth (wealth - substitute-capital)
          set owner? true
          ]
        ]
      ]
    ]



  let coal-patches patches with [pcolor = black]
  ask coal-patches [
    let coal-firms firms in-radius distance-setup
    ask coal-firms [
      set coal? true]
  ]


end
@#$#@#$#@
GRAPHICS-WINDOW
1114
10
1540
437
-1
-1
27.905
1
10
1
1
1
0
0
0
1
-7
7
-7
7
0
0
1
ticks
60.0

BUTTON
836
181
899
214
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
910
182
973
215
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
0
51
59
96
NIL
employed
2
1
11

MONITOR
67
54
139
99
NIL
unemployed
2
1
11

MONITOR
2
105
59
150
food prod.
aggregate-farm-production
2
1
11

MONITOR
0
156
60
201
good prod.
aggregate-industry-production
17
1
11

MONITOR
1
206
59
251
services
aggregate-services
17
1
11

MONITOR
66
104
123
149
metabolism
aggregate-metabolism
17
1
11

MONITOR
64
156
126
201
goods dem.
aggregate-goods-demand
2
1
11

MONITOR
65
206
126
251
serv. dem.
aggregate-services-demand
2
1
11

MONITOR
0
251
57
296
profits
aggregate-profits
17
1
11

MONITOR
58
251
123
296
dividends
aggregate-dividends
17
1
11

MONITOR
140
51
212
96
farms capital
average-capital-farms
2
1
11

MONITOR
214
50
296
95
firms capital
average-capital-firms
2
1
11

MONITOR
127
104
199
149
interest paid
yearly-interest
2
1
11

MONITOR
208
103
272
148
ind. firms
count firms with [industry? = true]
1
1
11

MONITOR
140
157
216
202
tot. fail firms
failed-firms
1
1
11

MONITOR
220
158
300
203
tot. new firms
new-firms
1
1
11

MONITOR
136
208
207
253
tot. bonds
bonds-sold
1
1
11

MONITOR
211
207
320
252
tot. new households
kids-number
1
1
11

MONITOR
135
255
206
300
avg. salary
mean-salaries
2
1
11

MONITOR
213
255
305
300
avg. farm price
mean-farm-price
2
1
11

MONITOR
133
304
213
349
avg. firm price
mean-firm-price
2
1
11

MONITOR
218
306
317
351
avg. service price
mean-service-price
2
1
11

MONITOR
0
305
57
350
NIL
taxes
2
1
11

MONITOR
59
303
127
348
NIL
GDP-income
2
1
11

MONITOR
0
353
78
398
NIL
GDP-spending
2
1
11

MONITOR
83
352
217
397
total goods sold value
goods-income-value
2
1
11

MONITOR
223
350
317
395
total value labor
labor-income-value
2
1
11

MONITOR
5
402
131
447
total food value sold
land-income-value
2
1
11

MONITOR
135
402
275
447
total value service sold
service-income-value
2
1
11

MONITOR
3
450
102
495
total value profits 
profit-income-value
2
1
11

MONITOR
109
450
228
495
NIL
mean-profits-farms
2
1
11

MONITOR
235
449
329
494
NIL
mean-profit-firms
2
1
11

TEXTBOX
38
19
253
49
Monitor for global macro variables
12
0.0
1

TEXTBOX
460
13
664
43
Choosers for Initial parametric values
12
0.0
1

TEXTBOX
338
60
488
78
number of breeds
11
0.0
1

TEXTBOX
327
319
477
337
industrial transition
11
0.0
1

TEXTBOX
352
177
502
195
initial prices
11
0.0
1

TEXTBOX
319
263
469
281
labor and productivities
11
0.0
1

TEXTBOX
366
462
516
480
others
11
0.0
1

MONITOR
0
675
77
720
NIL
labor-check
17
1
11

MONITOR
82
674
170
719
NIL
labor-check-1
17
1
11

MONITOR
176
674
264
719
NIL
labor-check-2
17
1
11

MONITOR
275
673
369
718
NIL
production-check
0
1
11

MONITOR
379
673
496
718
NIL
service-market-check
0
1
11

PLOT
0
721
582
1098
GDP
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"GDP income" 1.0 0 -14439633 true "plot GDP-income" "plot GDP-income"
"GDP spending" 1.0 0 -5298144 true "plot GDP-spending" "plot GDP-spending"
"real GDP" 1.0 0 -955883 true "" "plot real-GDP-spending"

PLOT
588
722
1188
1098
Shares of GDP
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Total income from goods" 1.0 0 -955883 true "plot goods-income-value" "plot goods-income-value"
"Total income from services" 1.0 0 -1184463 true "plot service-income-value" "plot service-income-value"
"Total income from land" 1.0 0 -13840069 true "plot land-income-value" "plot land-income-value"

PLOT
1191
722
1739
1099
Wealth
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"bourgeoisie" 1.0 0 -955883 true "plot mean [wealth] of bourgeoisie" "plot mean [wealth] of bourgeoisie"
"nobles" 1.0 0 -11221820 true "plot mean [wealth] of nobles" "plot mean [wealth] of nobles"
"workers" 1.0 0 -2674135 true "plot mean [wealth] of workers" "plot mean [wealth] of workers"

MONITOR
277
104
334
149
firms
count firms
0
1
11

MONITOR
0
499
116
544
average stock
mean [stock] of firms
2
1
11

MONITOR
124
499
242
544
potential industries
count firms with [stock = 0]
0
1
11

MONITOR
505
670
659
715
NIL
industrial-switching-check
17
1
11

MONITOR
667
671
845
716
NIL
mean [wealth] of government
17
1
11

TEXTBOX
358
385
508
403
government
11
0.0
1

CHOOSER
436
377
579
422
government-features
government-features
true false
0

CHOOSER
439
533
572
578
foreign-market-features
foreign-market-features
true false
1

MONITOR
851
671
978
716
NIL
government-check-1
17
1
11

MONITOR
981
671
1108
716
NIL
government-check-2
17
1
11

MONITOR
1112
672
1239
717
NIL
government-check-3
17
1
11

MONITOR
1247
671
1428
716
NIL
mean [needs] of government
17
1
11

MONITOR
27
621
145
666
NIL
firms-market-check
17
1
11

PLOT
0
1108
582
1479
Prices
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Service price" 1.0 0 -955883 true "plot mean [price] of bourgeoisie" "plot mean [price] of bourgeoisie"
"Goods price" 1.0 0 -13840069 true "plot mean [price] of firms" "plot mean [price] of firms"
"Food Price" 1.0 0 -14070903 true "plot mean [price] of farms" "plot mean [price] of farms"
"Labor Price" 1.0 0 -2674135 true "plot mean [labor-price] of workers" "plot mean [labor-price] of workers"

TEXTBOX
290
543
440
561
foreign-market (obsolete)
11
0.0
1

CHOOSER
579
533
709
578
initial-foreing-percentage
initial-foreing-percentage
0.05 0.1 0.15 0.2 0.3
0

CHOOSER
718
532
856
577
change-chance
change-chance
0.2 0.3 0.5 0.7
0

CHOOSER
862
533
1000
578
foreing-markup
foreing-markup
1 2 3 4
3

MONITOR
6
553
172
598
NIL
workers with [farm? = true]
17
1
11

INPUTBOX
438
41
543
101
number-of-workers
10000.0
1
0
Number

INPUTBOX
546
41
667
101
number-of-bourgeoisie
2000.0
1
0
Number

INPUTBOX
669
41
765
101
number-of-nobles
100.0
1
0
Number

INPUTBOX
772
41
861
101
number-of-firms
500.0
1
0
Number

INPUTBOX
868
41
959
101
number-of-cities
15.0
1
0
Number

INPUTBOX
965
40
1050
100
number-of-coal
5.0
1
0
Number

INPUTBOX
440
104
550
164
initial-labor-price
6.0
1
0
Number

TEXTBOX
342
119
492
137
initial endowments
11
0.0
1

INPUTBOX
555
105
681
165
initial-household-wealth
5.0
1
0
Number

INPUTBOX
686
105
811
165
bourgeoisie-ratio-wealth
30.0
1
0
Number

INPUTBOX
816
105
920
165
nobles-ratio-wealth
100.0
1
0
Number

INPUTBOX
925
104
1004
164
initial-capital
500.0
1
0
Number

INPUTBOX
437
169
543
229
initial-service-price
15.0
1
0
Number

INPUTBOX
548
170
608
230
markup
2.0
1
0
Number

INPUTBOX
614
171
729
231
price-change-chance
0.7
1
0
Number

INPUTBOX
734
171
802
231
price-delta
0.3
1
0
Number

INPUTBOX
437
237
536
297
firm-productivity
6.0
1
0
Number

INPUTBOX
543
237
653
297
service-productivity
10.0
1
0
Number

INPUTBOX
775
308
889
368
trans-cost-multiplier
250.0
1
0
Number

INPUTBOX
807
237
872
297
firm-labor
25.0
1
0
Number

INPUTBOX
878
238
982
298
firm-industrial-labor
100.0
1
0
Number

INPUTBOX
658
237
803
297
initial-industrial-productivity
15.0
1
0
Number

INPUTBOX
437
308
496
368
trans-cost
2500.0
1
0
Number

INPUTBOX
502
309
552
369
time
100.0
1
0
Number

INPUTBOX
561
309
620
369
safe-zone
0.3
1
0
Number

INPUTBOX
628
307
770
367
industrial-switch-probability
0.2
1
0
Number

INPUTBOX
583
377
696
437
government-demand
0.2
1
0
Number

INPUTBOX
702
376
846
436
percent-government-workers
0.05
1
0
Number

INPUTBOX
853
375
909
435
tax-rate
1.01
1
0
Number

INPUTBOX
915
374
1070
434
government-initial-wealth
10000.0
1
0
Number

INPUTBOX
437
446
506
506
distribution
7.0
1
0
Number

INPUTBOX
512
447
589
507
reproduction
250.0
1
0
Number

INPUTBOX
596
447
675
507
distance-setup
10.0
1
0
Number

INPUTBOX
684
446
819
506
service-distance-multiplier
1.5
1
0
Number

@#$#@#$#@
## MISSING

(Government demand for goods, food and services)
(foreign market is also missing



## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
